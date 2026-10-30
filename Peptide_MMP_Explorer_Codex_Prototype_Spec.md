# PepSeMan — スクラッチ実装指示書

更新: 2026-10-30

**公開境界と現行機能の追記:** `data/source/` の実CSVと `templates/sequence_dictionary.json` はローカル専用とし、GitHubに公開するソースツリーへ含めない。GitHub Pages公開は取りやめた。過去のGit履歴に残るデータの削除は別途承認を要する。現行実装は以下の直鎖中心の原型仕様を拡張し、head-to-tail、ジスルフィド、チオエーテル、環内C=Cに加えて、N置換1,2,3-トリアゾール架橋を構造から判定する。異なるtopology間をMMPにしない。細部は README と docs/ARCHITECTURE.md を参照する。

この文書は、空の作業ディレクトリから **現在の PepSeMan と同等のローカルアプリ**を実装するための受入仕様である。既存の SQLite DB、生成済みファイル、稼働中のサーバー、未記録の手作業に依存してはならない。ただし初期データの CSV 2 本は入力素材として必要であり、仕様書から生成してはならない。実装の背景と実測値は README.md、docs/ARCHITECTURE.md、docs/VALIDATION_REPORT.md を参照するが、矛盾する場合はこの文書の受入条件を優先し、差異を報告する。

## 1. 成果物と境界

- 直鎖ペプチドの SMILES/CSV/SMI/SDF を読み、5 残基以上で主鎖が一意に決まる分子だけを解析する。曖昧な化学構造は推測せず QC 理由付きで除外する。
- 残基または末端キャップ **1 箇所だけ**が異なる strict matched molecular pair (MMP) を抽出する。活性値がないデータでも構造ペアを扱う。
- MDM2 と GLP-1 のローカルCSVを別プロジェクトの初期データとして利用できるようにし、共通 Sequence dictionary を検証する。CSV自体は同梱しない。ランタイムで ChEMBL API に接続しない。
- UI は **Project / Analysis / Licenses の 3 タブ**とする。旧文書の独立した Import QC、R Decomposition、MMP Table、Pair Detail、Export、Job History タブは作らない。
- Web アプリは 127.0.0.1 のみに公開する。認証を実装していないため、外部インターフェースでの待受や公開プロキシへの接続は対象外。
- 分岐・複数フラグメントの自動解釈、残基挿入/削除、2 箇所以上の同時置換、活性予測、臨床判断は対象外。受理できない場合は「存在しない」と結論づけず unsupported と示す。対象の環化型は冒頭の追記と現行アーキテクチャを参照する。

## 2. 再現可能なリポジトリ

推奨構成は以下。モジュール分割は変えてよいが、入出力と受入条件は変えない。

~~~text
PepSeMan/
  app/                    FastAPI、化学解析、MMP、DB、worker
  frontend/               React/TypeScript/Vite UI、公開 SVG/PNG
  scripts/                初期データ投入、dictionary snapshot、macOS アイコン設定
  tests/                  化学、import、差分更新、API の pytest
  frontend/tests/         Playwright のブラウザテスト
  data/source/            ローカル専用の入力CSV（Git無視）
  templates/sequence_dictionary.json  ローカル専用（Git無視）
  assets/                 8 球の環状ペプチドを表すアイコン
  Launch PepSeMan.command
  PepSeMan.icns
  requirements.txt
  frontend/package-lock.json
  README.md
  docs/ARCHITECTURE.md
  docs/VALIDATION_REPORT.md
  docs/DEPENDENCIES.md
  LICENSE
  DATA_LICENSE.md
  .gitignore
~~~

Python 3.12、Node.js 22+、RDKit、FastAPI、Pydantic、SQLite WAL、React、TypeScript、TanStack Table、Plotly.js、pytest、Playwright を使う。Python 依存は requirements.txt、Node 依存は package-lock.json で固定する。現行の再現基準は RDKit 2026.9.1、FastAPI 0.136.1、Pydantic 2.13.4、Uvicorn 0.38.0、React 19.1.1、TanStack Table 8.21.3、Plotly.js 4.1.2、Vite 7.3.7、TypeScript 5.9.3、Playwright 1.55.1。残りの直接依存と推移依存も lock に固定し、更新が必要なら互換性とライセンスを再検証する。macOS Apple Silicon が既存の検証環境であり、他 OS は動作未検証と明記する。DuckDB 等の追加 DB は必要ない。UI ビルドは API が静的配信し、開発時のみ Vite の /api プロキシを使う。

**入力素材の確認:** MDM2、GLP-1、環状ペプチド、RAP1の実CSVは `data/source/` にローカル配置する。Git checkoutには含めない。必要なデータは権利を確認したうえで個別に入手し、欠落時はseedを失敗させる。既存DBやWebから代替データを推測生成して受入済みとしない。公開CIでは合成構造fixtureを主たる回帰テストに使い、実CSVの検証はファイルが存在するローカル環境で行う。

以下のコマンドがクリーンな checkout で実行できるようにする。初期投入は冪等とし、同一データの再投入で重複プロジェクトを作らない。

~~~sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cd frontend && npm ci && npm run build && cd ..
python -m app
~~~

別ターミナルで同じ仮想環境を有効化して python -m app.worker を起動する。ローカルCSVがある環境でのみ `python -m scripts.seed_prototype` と `python -m scripts.seed_ring_peptide` を実行する。標準の閲覧先は http://127.0.0.1:8000。README に起動/停止、環境変数 PMMP_DATA_DIR・PMMP_DB・PMMP_ASSET_DIR・PMMP_PORT、開発モード、テスト手順を記載する。実行DB、入力CSV、アップロード原本、ログ、PID、生成asset、dictionary snapshot、.venv、node_modules、dist、テスト成果物は .gitignore で除外する。

macOS の正式なランチャーは Launch PepSeMan.command とする。ローカル API と worker を起動しブラウザを開く。既存ポートの /api/health で同じデータディレクトリの instance と判定できる場合のみ再利用し、別 instance なら警告して停止する。ログ/PID は data/ に置く。Finder アイコン PepSeMan.icns、アプリ内アイコンと favicon は **8 個の球が円を作る環状ペプチドのモチーフ**に統一する。ランチャーは .venv を優先し、依存不足なら理由を表示する。ランチャー以外の重複起動スクリプトは増やさない。

## 3. テンプレートデータとライセンス

- MDM2.csv は ChEMBL 由来の CHEMBL5023 プロトタイプ。セミコロン区切り、49 列、492 測定行、404 ユニーク分子 ID。SMILES は 492/492 RDKit で読み取り可能という既存の基準を回帰テストに使う。
- GLP-1.csv は別プロジェクトの構造データ（現地ファイル 2131 データ行）。アッセイ/pChEMBL 列がないため、Value と ΔValue を捏造せず欠損で表示する。
- templates/sequence_dictionary.json は非天然残基の略称の **ローカル参照スナップショット**。Gitでは配布せず、ランタイムDBにも自動importしない。各構造から得たvariant_hashで共通辞書を構築し、構造式SVGは元分子から再生成する。
- LICENSE はアプリ本体・ドキュメント・自作アイコンを対象とする MIT（© 2026 keetane）。ChEMBL 由来 CSV と派生 dictionary は MIT の範囲外であり、DATA_LICENSE.md に出典、CC BY-SA 3.0、既知の由来不明点を明記する。元の ChEMBL release/download 詳細は記録されていないので推測で補わない。第三者 OSS は docs/DEPENDENCIES.md と Licenses タブで名称、用途、著作権者/出典、ライセンスとリンクを示す。

## 4. 入力、QC、化学構造

CSV、SMI、SDF を受け付ける。CSV はセミコロン区切りと quoted field に対応し、Smiles、Molecule ChEMBL ID、Standard Type/Relation/Value/Units、pChEMBL Value、Assay ChEMBL ID 等を列名で検出する。未知列を含む元行と原文 SMILES を保持する。1 分子に複数の測定行があるため、構造・プロジェクト所属・source record・activity を混同しない。SMI/SDF のプロパティとレコード ID も保持する。

RDKit で canonical isomeric SMILES を作り、標準化規則のバージョンと合わせて SHA-256 の structure_hash を計算する。元構造を勝手に中性化したり立体を補完しない。複数フラグメント/塩から主成分を推測で選ばない。SDF 往復は 2D 座標で十分であり、MDM2 の 492/492 レコードで canonical isomeric SMILES と非空プロパティの意味的同値性を検証する。3D 最適化は不要。

5mer 以上の直鎖ペプチドを、ペプチド結合の連なりから N→C 方向に **一意に**見つける。側鎖のアミドを主鎖に誤算入しない。残基ごとに N、Cα、carbonyl C/O、側鎖、主鎖 N の置換、立体、原子 index、前後の結合を記録する。R1/R2 は L/D に対応する Cα 側、R3 は主鎖 N 側の置換。Pro/D-Pro 等の環は R1/R2 と R3 に跨る bridge として一体で保持する。Gly、N-Me、ペプトイド、α,α-二置換、Hydroxyproline、側鎖修飾、N/C cap を保守的に判定する。L/D は R/S や画面上の向きから機械的に決めず、確定できなければ unknown とする。cap は残基数に含めない。

QC status は accepted、excluded_short_peptide、excluded_non_peptide、excluded_non_linear、excluded_multiple_fragments、unsupported_backbone 等を区別し、具体的な exclusion_reason を保存する。1–4mer は短鎖として除外、曖昧/分岐/環状は勝手な主鎖を採用しない。解析で使えなかったソース行も QC 一覧から消さない。

## 5. 配列表記と共有 dictionary

天然 L アミノ酸は大文字一文字、天然 D アミノ酸は小文字一文字、N 末 Ac は [Ac] とする。D-Ile/D-Thr のような側鎖キラリティも区別し、Gly に D 略号は作らない。これらの固定コードは編集対象から除外する。その他の非天然 building block は [4FP] 等の bracketed 略称を variant_hash ごとに一意に割り当て、名称未設定でも決定的な hash 付き略号にフォールバックする。衝突、括弧の混入、空文字、過長名を検証する。C cap などを表示する際も別残基として数えない。

編集可能なのは **bracketed 非天然残基だけ**。global sequence_aliases を variant_hash で管理し、変更した略称は MDM2 と GLP-1 の両方の R 分解、MMP table、Pair detail、CSV export に即時反映する。変更のために分解/MMP を再計算しない。編集値は括弧を除く 1–16 文字の ASCII 英数字・_・-、先頭は英数字とし、天然一文字コード、Ac/H/OH/NH2、既存の自動/手動略称との衝突を拒否する。空値は削除ではなく default 略称へのリセットとして扱う。辞書一覧には略称、default/custom 値、残基タイプ、立体、R1/R2/R3/bridge、元分子から描いた当該 building block の構造式、検索、保存/リセットを備える。snapshot の更新は明示的な開発コマンドで行い、通常の UI 編集で追跡ファイルを書き換えない。

## 6. strict MMP と活性差

同じ残基数の分子について、N cap・順序付き残基・C cap の各位置を 1 つだけマスクした context_hash を作る。同じ project・parser/engine version の同一 context bucket 内だけで variant_hash が違う組を検証する。変更部位以外は立体込みで一致、変更部位は残基または cap の 1 箇所のみ。挿入/削除や 2 箇所同時変更は不採用。全分子 O(N²) の総当たりを基本方式にしない。

pair_key は構造 hash を安定順に並べて作り、逆向きの重複保存を防ぐ。表示時の選択 parent/reference に合わせて Before/After、Original/Convert、ΔMW、ΔValue の符号を反転する。structure MMP と測定結合は分離し、測定のないペアも残す。

Standard Relation の引用符付き値を正規化しつつ原文も保持する。数値の ΔpChEMBL は **双方 relation =、pChEMBL が有効、同一 assay ID、同一 endpoint** のときだけ after − before とする。< や >、異なる測定条件、欠損値に厳密な数値差を与えない。IC50/Ki/Kd を混ぜない。GLP-1 の活性値は欠損表示とする。分子量は構造から算出し、ΔMW は同じ表示方向で after − before。

## 7. SQLite、ジョブ、差分更新

SQLite は WAL・foreign keys・明示 transaction・schema_migrations を使う。少なくとも projects、imports、source_records、molecules、project_molecules、activities、decompositions、residues、sequence_aliases、contexts、pairs、jobs、analysis_revisions を分ける。depiction cache は描画設定と構造 version で無効化できる設計とする。decomposition は structure_hash + parser_version + standardization_version で全プロジェクト共用、context/pair/測定/履歴は project 単位。source の原本と provenance は revision 復元後も追跡できるようにする。

import は API が原本をローカル保存して durable job を queue し、独立 worker が処理する。ブラウザを閉じても動き、再起動時にキュー/異常終了を安全に扱う。解析は DB 書込 transaction の外で計算し、完了した revision の有効化とソース/測定/context/pair の反映を原子的に commit する。失敗時は前の revision を維持し、エラーを job に記録する。復元は履歴を巻き戻さず、新しい revision として queue する。

同一ファイル SHA-256 と解析 version の再 import は重複として skip。未変更構造は分解を再計算せず、追加・修正・削除が影響する context bucket だけ再評価する。測定値だけの変更なら構造解析と MMP 再計算は 0、描画だけの変更なら depiction cache のみ無効化する。バージョン変更で安全な依存範囲を特定できなければ、その project を再解析する。activity ID がない測定行を不確かな同一性で上書きしない。ジョブには added/removed/changed/recomputed/reused と経過時間を残す。

project 名は作成日の Asia/Tokyo 日付を付けた **yyyymmdd_filename** とし、拡張子を除いた basename を安全に正規化し、重複時は連番 suffix を付ける。既存 DB があれば migration で既存名もこの規則に揃え、project ID とデータは保持する。

## 8. API 契約

FastAPI は /api/health を返し、全 project ID、reference_hash、pair_key を所属 project と照合する。他 project の URL を指定しても内容を返さず 404 とする。page は 1 以上、page_size は最大 500、sort は許可済み列だけを受け、値は SQL に直接埋め込まない。sort は表示列と同じ方向・型で全件に適用してから page を切り、NULL は昇順・降順とも末尾、同順位は pair_key 等の安定キー順とする。不正な page/sort/入力形式には 422 を返す。表の検索入力、略称、ファイル名、CSV/SMILES/SDF は境界で検証する。

最低限の route と用途:

| Route | 用途 |
| --- | --- |
| GET/POST /api/projects; PATCH/DELETE /api/projects/{id} | project 管理 |
| POST/GET /api/projects/{id}/imports | 原本投入・履歴 |
| GET /api/projects/{id}/qc | status/ID 検索、sort、page |
| GET /api/projects/{id}/decompositions | value_metric、sequence_filter、all_rows、sort |
| GET /api/projects/{id}/pairs | reference_hash、sort、page、配列と変更 span |
| GET /api/projects/{id}/pairs/selection | 表示ページを跨ぐ全選択用 pair key |
| GET /api/projects/{id}/pairs/{pair_key} | 方向付き pair detail |
| GET /api/projects/{id}/pairs/{pair_key}/preview | 小型 SVG |
| POST /api/projects/{id}/export | parent + 選択 partner の CSV |
| GET /api/projects/{id}/jobs、/revisions; POST /revisions/{revision_id}/restore | 永続 job / revision |
| GET /api/sequence-dictionary; PUT /api/sequence-dictionary/{variant_hash}; GET /api/sequence-dictionary/{variant_hash}/svg | 共通略称一覧・編集・構造式 |

一覧応答は rows、total、page、page_size を含む。R 分解は ID、nPairs、Value、sequence、構造・残基列を返す。Value は pChEMBL Value を既定とし、分子量/残基数へ切替可能。sequence_filter は配列表記を対象とする大文字小文字区別検索。MMP 応答は parent_sequence、parent_sequence_highlight_start/end、original、convert、site_position、delta_activity、delta_mw、pair_key、両側の構造 hash、活性/assay/provenance を返す。highlight は表示配列の **先頭 0・終端排他的な文字オフセット**であり、bracketed トークン全体を囲む。/pairs/selection は page に依存しない pair_key、first_structure_hash、second_structure_hash の配列を返す。pair detail は query の reference_hash で方向を指定し、parent の scaffold SVG と変更部位だけの Original/Convert BB SVG、before/after と property 差を返す。

CSV export の body は mode = selected_compounds と structure_hashes = [現在の parent, 選択 pair の partner 構造 ...] とする。UI から pair key 自体を export body に渡さず、partner hash に解決する。API は各 hash が当該 project に属することを再検証し、指定外の構造や他 project の source/測定を出さない。parent を必ず含め、source ID、canonical SMILES、sequence、MW、測定値、relation、assay、provenance を失わない。出力は **source ID と測定行の粒度**を保持し、同じ構造に複数の source ID/測定があれば必要な行数で表す。構造 hash の重複除去によって測定を落とさない。測定が無い構造も少なくとも 1 行は出す。

## 9. 画面と操作

タブバー左端に 8 球アイコン + PepSeMan、中央に Project / Analysis / Licenses、右端に Active project セレクタを置く。別の大きなタイトル header は作らない。URL の project/reference/pair を読み書きし、存在しない組合せでは安全に初期選択へ戻す。

### Project

Project 作成/選択、Import & QC、Job History / revision、Sequence dictionary を **同じタブ**に並べる。QC records は 1 ページ 5 行、ID 部分一致検索、page jump、sort を備える。取り込み進捗/除外理由/再 import の結果を確認できる。辞書は project を選び直しても同じ global alias を編集する。

### Analysis

- 左約 1/4 に R Decomposition。独立した説明文は不要。固定列は **ID、nPairs、Value** の順で左寄せし、残りの列を横スクロールする。section header 右端の Value セレクタは既定 pChEMBL Value。sequence 検索フィルタを置き、page 分割を行わない。
- 右約 3/4 は縦 2 分割。**上に Pair detail、下に MMP Table / Scatter Plot の切替 section**を置く。詳細には Parent sequence と parent scaffold、Original BB、Convert BB、ΔChEMBL Value/ΔMW 等を表示する。初期表示では R 分解の先頭 accepted 行、その parent の MMP table の先頭行を選ぶ。両テーブルは行クリックとフォーカス時の上下キーで行を切り替え、選択行を可視範囲に保つ。選択した reference と pair は詳細/散布図に同期する。
- Parent scaffold は parent 分子の N 末を左、C 末を右にした RDKit SVG。主鎖は自然なジグザグを保ちつつ概ね横方向、側鎖は可能な範囲で上下交互。proline 環等の制約を壊してまで交互にしない。分子グラフや stereochemical wedge を鏡像反転・改変しない。変換部位を黄色で highlight。整列不能なら警告と通常構造図を示す。
- scaffold 上の「N末→C末」見出しは置かず、元の **Parent sequence** を表示する。変換された residue/cap のトークンだけ黄色にする。MMP table の Parent sequence も -X- に置換せず同じ強調表示。Original/Convert は変換があった building block の短い配列表記とし、BB SVG も当該アミノ酸部分のみ。cap 差分なら cap の情報を明示する。
- MMP Table の先頭は Select checkbox。ヘッダー checkbox は現在の parent の全 pair をページ跨ぎで選択/解除。以降の列は **Parent sequence、Original、Convert、Site position、ΔValue、ΔMW** の順を先頭に置き、その後 ID、活性、assay、構造情報等を続ける。列幅をドラッグ/キーボードで変更可能とし、Parent sequence は既定で全文を読める幅にする。すべての実データ列を型に沿って sort 可能とし、選択/ページ状態を壊さない。MMP table 自体の大きな header は不要。Export parent + selected molecules ボタンは **MMP/Scatter のタブバー内**に置き、別 Export タブは作らない。
- Scatter は MW/ΔMW/活性/Δ活性から X/Y を選択でき、両軸の数値がある pair のみ点を描く。点クリックで pair detail を選ぶ。GLP-1 のような活性欠損データでも MW 軸で使える。

### Licenses

コード、同梱データ、主要 OSS のクレジット/用途/ライセンス/出典リンクを表で表示し、LICENSE と DATA_LICENSE.md の適用範囲を取り違えない。

## 10. 必須テストと実測受入

pytest で少なくとも以下を検証する。

1. 化学: 4mer 除外/5mer 採用、Ac/NH2 cap、L↔D、Ala↔Val、Pro/D-Pro/Hyp、N-Me、ペプトイド、Gly、修飾 Lys/Ser、α,α-二置換、側鎖アミド、未指定立体、複数フラグメント、環状/分岐を区別。同一構造の別 SMILES は同じ hash。2 箇所差は MMP にしない。
2. import/QC: MDM2 全 492 行を保持し、原本/測定 492 とユニーク分子 ID 404 を混同しない。CSV→2D SDF→readback の 492/492 構造・非空 property の意味的同値性。CSV/SDF 経由の受理/分解/MMP 結果一致。GLP-1 の活性欠損を確認。
3. 差分/永続化: 同一ファイル再投入は重複を作らず、活性だけ変更で構造再解析 0、追加/削除/構造変更は影響 bucket のみ更新、完全再計算と結果一致。途中失敗で旧 revision を維持し、復元は新 revision。worker を API/ブラウザと独立させる。2 project 間で source/測定/pair を漏らさず、dictionary alias だけ共用。
4. API/UI: project 外 reference/pair 拒否、sort と pagination の一貫性、QC の ID 検索/5 行/page jump、R 分解無 pagination・sequence filter・Value 切替、両表の上下キー、MMP sort/列幅/全選択/CSV、highlight の位置、散布図から詳細への選択、3 タブ、Licenses リンクを検証。

既存 MDM2 の参考実測値は **accepted 50 source rows / 49 unique 5mer+ structures、29 strict MMP pairs、同一 assay/endpoint の厳密値比較が可能な pair 9**。これは 2026-10-28 の parser とデータでの測定値であり、旧ヒューリスティックの「68 行」は正式な受理基準ではない。実装後に受理/除外内訳、GLP-1 の件数、処理時間と cache 再利用を実測し docs/VALIDATION_REPORT.md に記録する。数値が違えばテストを緩めず原因と仕様差を報告する。

受入コマンド:

~~~sh
python -m pytest
cd frontend && npm run build && npm run test:e2e
~~~

Playwright のブラウザ実体が未導入なら、その点を未実施と明記する。コード変更に応じて lint/型チェック、依存の脆弱性・ライセンス確認を行う。テストに通っていない機能を完成と報告しない。

## 11. 実装順と報告

1. 空ディレクトリ、依存 lock、設定、DB schema/migration、fixture、README の最小起動経路を作る。
2. import/QC、保存、SDF 往復、保守的な N→C 分解、化学ユニットテストを完成させる。
3. structure cache、context index、strict MMP、activity 結合、差分更新、worker/revision と API テストを完成させる。
4. global dictionary、3 タブ UI、描画、sort/keyboard/export、Playwright を完成させる。
5. 2 つの同梱 CSV をクリーンな DB に seed し、実測値/制約/ライセンスを文書化する。ランチャー/アイコン、.gitignore、配布前の秘密情報混入確認を行う。

各段階で「変更点、実際に実行した確認と結果、未確認事項、安全上の懸念」を簡潔に報告する。既存の未コミット変更と個人データを保護し、必要な最小差分だけ編集する。**git commit、push、公開、release/deploy はこの文書だけを根拠に自動実行しない。** それらは直近のユーザー依頼または明示の許可がある場合に、対象と影響を確認して行う。
