#!/bin/zsh
set -u

PROJECT_DIR="${0:A:h}"
if [[ -f "$PROJECT_DIR/PepSeMan.icns" && -x /usr/bin/swift ]] && ! /usr/bin/xattr -p com.apple.ResourceFork "$0" >/dev/null 2>&1; then
  /usr/bin/swift "$PROJECT_DIR/scripts/set_launcher_icon.swift" "$PROJECT_DIR/PepSeMan.icns" "${0:A}" >/dev/null 2>&1 || true
fi
DATA_DIR="$PROJECT_DIR/data"
PORT="${PMMP_PORT:-48001}"
APP_URL="http://127.0.0.1:$PORT/?project=1&tab=Project"
LOG_FILE="$DATA_DIR/launcher.log"
SERVER_PID_FILE="$DATA_DIR/server.pid"
WORKER_PID_FILE="$DATA_DIR/worker.pid"

mkdir -p "$DATA_DIR"

find_python() {
  local candidate
  for candidate in \
    "$PROJECT_DIR/.venv/bin/python" \
    "/Users/keetane/miniforge3/bin/python" \
    "/opt/homebrew/bin/python3" \
    "/usr/local/bin/python3" \
    "/usr/bin/python3"; do
    if [[ -x "$candidate" ]] && "$candidate" -c 'import fastapi, rdkit, uvicorn' >/dev/null 2>&1; then
      print -r -- "$candidate"
      return 0
    fi
  done
  return 1
}

pid_is_running() {
  local pid_file="$1"
  [[ -f "$pid_file" ]] || return 1
  local pid
  pid="$(<"$pid_file")"
  [[ "$pid" == <-> ]] && kill -0 "$pid" >/dev/null 2>&1 || return 1
  [[ "$(/bin/ps -p "$pid" -o command= 2>/dev/null)" == *' -m app.worker'* ]]
}

server_is_ready() {
  local response
  response="$(/usr/bin/curl --silent --fail --max-time 1 "http://127.0.0.1:$PORT/api/health" 2>/dev/null)" || return 1
  [[ "$response" == *'"status":"ok"'* && "$response" == *'"instance_id":"'"$EXPECTED_INSTANCE"'"'* ]]
}

show_error() {
  print -u2 -- "$1"
  /usr/bin/osascript -e 'on run argv' -e 'display dialog (item 1 of argv) with title "PepSeMan" buttons {"OK"} default button "OK" with icon stop' -e 'end run' "$1" >/dev/null 2>&1
}

PYTHON_BIN="$(find_python)" || {
  show_error "Python環境が見つかりません。READMEの手順で依存関係をインストールしてください。"
  exit 1
}

cd "$PROJECT_DIR" || exit 1
EXPECTED_INSTANCE="$(PMMP_DATA_DIR="$DATA_DIR" "$PYTHON_BIN" -c 'import hashlib; from app.config import settings; print(hashlib.sha256(str(settings.database_path).encode()).hexdigest()[:16])')" || exit 1

if /usr/bin/curl --silent --fail --max-time 1 "http://127.0.0.1:$PORT/api/health" >/dev/null 2>&1 && ! server_is_ready; then
  show_error "Port $PORT は別のデータを使うサーバーが利用中です。古いサーバーを終了してから再起動してください。"
  exit 1
fi

if ! server_is_ready; then
  PMMP_PORT="$PORT" PMMP_DATA_DIR="$DATA_DIR" "$PYTHON_BIN" -m app </dev/null >>"$LOG_FILE" 2>&1 &
  STARTED_SERVER_PID="$!"
  print -r -- "$STARTED_SERVER_PID" > "$SERVER_PID_FILE"
  integer attempt=0
  while (( attempt < 40 )) && ! server_is_ready; do
    /bin/sleep 0.25
    (( attempt++ ))
  done
fi

if ! server_is_ready; then
  show_error "アプリサーバーを起動できませんでした。ログ: $LOG_FILE"
  exit 1
fi

if ! pid_is_running "$WORKER_PID_FILE"; then
  PMMP_DATA_DIR="$DATA_DIR" "$PYTHON_BIN" -m app.worker </dev/null >>"$LOG_FILE" 2>&1 &
  print -r -- "$!" > "$WORKER_PID_FILE"
fi

if [[ "${PMMP_SKIP_OPEN:-0}" != "1" ]]; then
  /usr/bin/open "$APP_URL"
fi

if [[ -n "${STARTED_SERVER_PID:-}" ]]; then
  wait "$STARTED_SERVER_PID"
fi
