import { defineConfig } from "@playwright/test";
import { mkdtempSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";

const dataDir = process.env.PMMP_E2E_DATA_DIR ?? mkdtempSync(join(tmpdir(), "pepseman-e2e-"));

export default defineConfig({
  testDir: "./tests",
  use: { baseURL: "http://127.0.0.1:48000" },
  webServer: {
    command: "python -m app",
    url: "http://127.0.0.1:48000/api/health",
    cwd: "..",
    reuseExistingServer: false,
    env: { PMMP_PORT: "48000", PMMP_DATA_DIR: dataDir },
  },
});

