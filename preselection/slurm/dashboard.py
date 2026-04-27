#!/usr/bin/env python3
"""
Live web dashboard for SLURM job status.

Serves an auto-refreshing HTML page showing per-task status with progress bars.

Usage (run from preselection/):
    python3 slurm/dashboard.py --date 20260416
    python3 slurm/dashboard.py --date 20260416 --port 8080
    python3 slurm/dashboard.py --date 20260416 --refresh 15

Since this runs on a remote login node, you need an SSH tunnel to view the
dashboard in your local browser. Open a separate terminal and run:

    ssh -L 8085:localhost:8085 <user>@<host>

Then open http://localhost:8085 in your browser.
"""

import argparse
import json
import os
import sys
from collections import defaultdict
from http.server import HTTPServer, BaseHTTPRequestHandler
from pathlib import Path
from urllib.parse import urlparse, parse_qs

# Import status.py functions from the same directory
sys.path.insert(0, str(Path(__file__).parent))
from status import (
    load_manifest, get_slurm_job_states, determine_job_status, SLURM_OUTPUT_DIR
)

SCRIPT_DIR = Path(__file__).parent.resolve()
JOBS_DIR = SCRIPT_DIR / SLURM_OUTPUT_DIR


def get_tasks_for_date(date_str: str) -> list[str]:
    """Find all task directories matching a date string."""
    if not JOBS_DIR.exists():
        return []
    return sorted([
        d.name for d in JOBS_DIR.iterdir()
        if d.is_dir() and date_str in d.name
    ])


def get_all_status(date_str: str) -> dict:
    """Query SLURM and return status for all tasks matching a date."""
    tasks = get_tasks_for_date(date_str)
    results = []

    grand = {"total": 0, "completed": 0, "running": 0, "queued": 0, "failed": 0}

    for task_name in tasks:
        task_dir = JOBS_DIR / task_name
        manifest = load_manifest(task_dir)
        if not manifest:
            continue

        jobs = manifest.get("jobs", {})
        if not jobs:
            continue

        # Query SLURM
        slurm_ids = [
            j["slurm_job_id"] for j in jobs.values()
            if j.get("slurm_job_id") is not None
        ]
        slurm_states, _ = get_slurm_job_states(slurm_ids)

        # Count statuses
        counts = defaultdict(int)
        for job_id, job in jobs.items():
            status = determine_job_status(job, slurm_states)
            counts[status] += 1

        n_total = len(jobs)
        n_completed = counts.get("completed", 0)
        n_running = counts.get("running", 0)
        n_queued = counts.get("queued", 0) + counts.get("submitted", 0)
        n_failed = counts.get("failed", 0)

        # Strip prefix and date suffix for short name
        short_name = task_name
        if short_name.startswith("merged_"):
            short_name = short_name[len("merged_"):]
        # Remove _<date>_<channel> suffix
        idx = short_name.find(f"_{date_str}")
        if idx >= 0:
            short_name = short_name[:idx]

        results.append({
            "task": task_name,
            "short_name": short_name,
            "total": n_total,
            "completed": n_completed,
            "running": n_running,
            "queued": n_queued,
            "failed": n_failed,
        })

        grand["total"] += n_total
        grand["completed"] += n_completed
        grand["running"] += n_running
        grand["queued"] += n_queued
        grand["failed"] += n_failed

    return {"tasks": results, "grand": grand, "n_tasks": len(results)}


def make_html(date_str: str, refresh_sec: int) -> str:
    """Generate the dashboard HTML page."""
    return f"""\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SLURM Dashboard - {date_str}</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{
    background: #1a1a2e; color: #e0e0e0;
    font-family: 'SF Mono', 'Fira Code', 'Consolas', monospace;
    padding: 24px;
  }}
  h1 {{ color: #e94560; margin-bottom: 4px; font-size: 1.4em; }}
  .subtitle {{ color: #666; margin-bottom: 20px; font-size: 0.85em; }}
  #status-msg {{
    color: #888; font-size: 0.8em; margin-bottom: 12px;
  }}
  table {{
    width: 100%; border-collapse: collapse;
    font-size: 0.85em;
  }}
  th {{
    text-align: left; padding: 8px 12px;
    border-bottom: 2px solid #333;
    color: #888; font-weight: 600;
    position: sticky; top: 0; background: #1a1a2e;
  }}
  th.num {{ text-align: right; }}
  td {{
    padding: 6px 12px; border-bottom: 1px solid #222;
  }}
  td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  tr:hover {{ background: #16213e; }}
  tr.grand {{ border-top: 2px solid #444; font-weight: bold; }}
  tr.grand td {{ padding-top: 12px; }}

  .done {{ color: #4ecca3; }}
  .run  {{ color: #4ea8de; }}
  .queue {{ color: #f0c040; }}
  .fail {{ color: #e94560; }}
  .zero {{ color: #444; }}

  .bar-container {{
    width: 160px; height: 16px;
    background: #222; border-radius: 3px;
    overflow: hidden; display: inline-block;
    vertical-align: middle;
  }}
  .bar-fill {{
    height: 100%; border-radius: 3px;
    transition: width 0.5s ease;
  }}
  .bar-fill.green {{ background: linear-gradient(90deg, #2d6a4f, #4ecca3); }}
  .bar-fill.blue  {{ background: linear-gradient(90deg, #1a4a6e, #4ea8de); }}
  .bar-fill.red   {{ background: linear-gradient(90deg, #6e1a1a, #e94560); }}

  .pct {{
    display: inline-block; width: 42px;
    text-align: right; margin-left: 8px;
    font-size: 0.85em; color: #888;
  }}

  .row-done td:first-child {{ border-left: 3px solid #4ecca3; }}
  .row-run  td:first-child {{ border-left: 3px solid #4ea8de; }}
  .row-fail td:first-child {{ border-left: 3px solid #e94560; }}
  .row-queue td:first-child {{ border-left: 3px solid #f0c040; }}

  .spinner {{
    display: inline-block; width: 12px; height: 12px;
    border: 2px solid #444; border-top-color: #4ea8de;
    border-radius: 50%; animation: spin 0.8s linear infinite;
    margin-right: 6px; vertical-align: middle;
  }}
  @keyframes spin {{ to {{ transform: rotate(360deg); }} }}
</style>
</head>
<body>

<h1>SLURM Job Dashboard</h1>
<div class="subtitle">Date filter: {date_str} &mdash; Refresh: {refresh_sec}s</div>
<div id="status-msg">Loading...</div>

<table>
<thead>
  <tr>
    <th>Task</th>
    <th class="num">Total</th>
    <th class="num">Done</th>
    <th class="num">Run</th>
    <th class="num">Queue</th>
    <th class="num">Fail</th>
    <th>Progress</th>
  </tr>
</thead>
<tbody id="tbody"></tbody>
</table>

<script>
const REFRESH = {refresh_sec} * 1000;
const DATE = "{date_str}";

function cls(n, c) {{ return n > 0 ? c : 'zero'; }}

function rowClass(t) {{
  if (t.total > 0 && t.completed === t.total) return 'row-done';
  if (t.failed > 0 && t.running === 0 && t.queued === 0) return 'row-fail';
  if (t.running > 0 || t.queued > 0) return 'row-run';
  return 'row-queue';
}}

function barClass(t) {{
  if (t.total > 0 && t.completed === t.total) return 'green';
  if (t.failed > 0 && t.running === 0 && t.queued === 0) return 'red';
  return 'blue';
}}

function render(data) {{
  const tbody = document.getElementById('tbody');
  let html = '';

  for (const t of data.tasks) {{
    const pct = t.total > 0 ? Math.round(100 * t.completed / t.total) : 0;
    const rc = rowClass(t);
    const bc = barClass(t);

    html += `<tr class="${{rc}}">
      <td>${{t.short_name}}</td>
      <td class="num">${{t.total}}</td>
      <td class="num ${{cls(t.completed,'done')}}">${{t.completed}}</td>
      <td class="num ${{cls(t.running,'run')}}">${{t.running}}</td>
      <td class="num ${{cls(t.queued,'queue')}}">${{t.queued}}</td>
      <td class="num ${{cls(t.failed,'fail')}}">${{t.failed}}</td>
      <td>
        <div class="bar-container">
          <div class="bar-fill ${{bc}}" style="width:${{pct}}%"></div>
        </div>
        <span class="pct">${{pct}}%</span>
      </td>
    </tr>`;
  }}

  // Grand total row
  const g = data.grand;
  const gpct = g.total > 0 ? Math.round(100 * g.completed / g.total) : 0;
  html += `<tr class="grand">
    <td>TOTAL (${{data.n_tasks}} tasks)</td>
    <td class="num">${{g.total}}</td>
    <td class="num done">${{g.completed}}</td>
    <td class="num run">${{g.running}}</td>
    <td class="num queue">${{g.queued}}</td>
    <td class="num fail">${{g.failed}}</td>
    <td>
      <div class="bar-container">
        <div class="bar-fill blue" style="width:${{gpct}}%"></div>
      </div>
      <span class="pct">${{gpct}}%</span>
    </td>
  </tr>`;

  tbody.innerHTML = html;
}}

async function refresh() {{
  const msg = document.getElementById('status-msg');
  msg.innerHTML = '<span class="spinner"></span>Querying SLURM...';
  try {{
    const resp = await fetch(`/api/status?date=${{DATE}}`);
    const data = await resp.json();
    render(data);
    const now = new Date().toLocaleTimeString();
    msg.textContent = `Last updated: ${{now}} — Next refresh in ${{REFRESH/1000}}s`;
  }} catch (e) {{
    msg.textContent = `Error: ${{e.message}}`;
  }}
}}

refresh();
setInterval(refresh, REFRESH);
</script>
</body>
</html>"""


class DashboardHandler(BaseHTTPRequestHandler):
    """HTTP request handler for the dashboard."""

    date_str = ""
    refresh_sec = 30

    def log_message(self, format, *args):
        # Quieter logging
        sys.stderr.write(f"[{self.date_time_string()}] {format % args}\n")

    def do_GET(self):
        parsed = urlparse(self.path)

        if parsed.path == "/":
            html = make_html(self.date_str, self.refresh_sec)
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(html.encode())

        elif parsed.path == "/api/status":
            params = parse_qs(parsed.query)
            date = params.get("date", [self.date_str])[0]
            data = get_all_status(date)
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Cache-Control", "no-cache")
            self.end_headers()
            self.wfile.write(json.dumps(data).encode())

        else:
            self.send_response(404)
            self.end_headers()


def main():
    parser = argparse.ArgumentParser(description="SLURM job status dashboard")
    parser.add_argument("--date", "-d", required=True,
                        help="Date string to filter tasks (e.g. 20260416)")
    parser.add_argument("--port", "-p", type=int, default=8085,
                        help="Port to serve on (default: 8085)")
    parser.add_argument("--refresh", "-r", type=int, default=30,
                        help="Auto-refresh interval in seconds (default: 30)")
    args = parser.parse_args()

    DashboardHandler.date_str = args.date
    DashboardHandler.refresh_sec = args.refresh

    # Verify tasks exist
    tasks = get_tasks_for_date(args.date)
    if not tasks:
        print(f"No tasks found matching date '{args.date}'")
        sys.exit(1)

    print(f"Found {len(tasks)} tasks matching '{args.date}'")
    print(f"Dashboard: http://localhost:{args.port}")
    print(f"Refresh interval: {args.refresh}s")
    print("Press Ctrl+C to stop.\n")

    server = HTTPServer(("0.0.0.0", args.port), DashboardHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == "__main__":
    main()
