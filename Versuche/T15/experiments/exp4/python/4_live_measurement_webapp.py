import csv
import threading
import time
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime
from typing import Deque, Dict, List, Optional

import board
import busio
import adafruit_lsm9ds1
import RPi.GPIO as GPIO
from flask import Flask, jsonify, Response


# ============================================================
# Konfiguration
# ============================================================

HOST = "0.0.0.0"
PORT = 5000

SWITCH_GPIO = 14
CSV_FILENAME = "live_data_new.csv"

# Browser aktualisiert diese API regelmäßig
API_REFRESH_INTERVAL_MS = 500

# Messintervall
MEASUREMENT_INTERVAL_SECONDS = 0.05

# Anzahl Punkte, die im RAM gehalten werden
MAX_POINTS = 3000

# Anzahl Zeilen in der Tabelle
TABLE_POINTS = 10

# Zeitfenster des Plots in Sekunden
PLOT_WINDOW_SECONDS = 45.0


# ============================================================
# Datenmodell
# ============================================================

@dataclass
class MeasurementPoint:
    timestamp_unix: float
    timestamp_iso: str
    accel_x: float
    accel_y: float
    accel_z: float
    gyro_x: float
    gyro_y: float
    gyro_z: float
    mag_x: float
    mag_y: float
    mag_z: float


@dataclass
class AppState:
    is_measuring: bool = False
    status_text: str = "Initialisiert"
    last_error: Optional[str] = None
    last_data_timestamp_unix: Optional[float] = None
    points: Deque[MeasurementPoint] = field(default_factory=lambda: deque(maxlen=MAX_POINTS))
    lock: threading.Lock = field(default_factory=threading.Lock)

    def add_point(self, point: MeasurementPoint) -> None:
        with self.lock:
            self.points.append(point)
            self.last_data_timestamp_unix = point.timestamp_unix
            self.is_measuring = True
            self.status_text = "Misst"

    def set_idle(self, text: str = "Nicht am Messen") -> None:
        with self.lock:
            self.is_measuring = False
            self.status_text = text

    def set_error(self, message: str) -> None:
        with self.lock:
            self.is_measuring = False
            self.status_text = "Fehler"
            self.last_error = message

    def clear_error(self) -> None:
        with self.lock:
            self.last_error = None

    def snapshot(self) -> Dict:
        with self.lock:
            points_list = list(self.points)

            return {
                "is_measuring": self.is_measuring,
                "status_text": self.status_text,
                "last_error": self.last_error,
                "last_data_timestamp_unix": self.last_data_timestamp_unix,
                "table_points": [
                    {
                        "timestamp_unix": p.timestamp_unix,
                        "timestamp_iso": p.timestamp_iso,
                        "accel_x": p.accel_x,
                        "accel_y": p.accel_y,
                        "accel_z": p.accel_z,
                        "gyro_x": p.gyro_x,
                        "gyro_y": p.gyro_y,
                        "gyro_z": p.gyro_z,
                        "mag_x": p.mag_x,
                        "mag_y": p.mag_y,
                        "mag_z": p.mag_z,
                    }
                    for p in points_list[-TABLE_POINTS:]
                ],
                "plot_points": [
                    {
                        "timestamp_unix": p.timestamp_unix,
                        "timestamp_iso": p.timestamp_iso,
                        "accel_z": p.accel_z,
                    }
                    for p in points_list
                ],
                "point_count": len(points_list),
            }


# ============================================================
# Globale Objekte
# ============================================================

APP = Flask(__name__)
STATE = AppState()

sensor = None
csv_file = None
csv_writer = None


# ============================================================
# Hilfsfunktionen
# ============================================================

def format_exception_message(exc: Exception) -> str:
    return f"{exc.__class__.__name__}: {exc}"


def unix_to_iso(ts: float) -> str:
    return datetime.fromtimestamp(ts).isoformat(timespec="milliseconds")


# ============================================================
# Hardware / CSV Setup
# ============================================================

def initialize_hardware() -> None:
    global sensor

    GPIO.setmode(GPIO.BCM)
    GPIO.setup(SWITCH_GPIO, GPIO.IN, pull_up_down=GPIO.PUD_UP)

    i2c = busio.I2C(board.SCL, board.SDA)
    sensor = adafruit_lsm9ds1.LSM9DS1_I2C(i2c)


def initialize_csv() -> None:
    global csv_file, csv_writer

    csv_file = open(CSV_FILENAME, "w", newline="")
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow([
        "Time",
        "accel_x", "accel_y", "accel_z",
        "gyro_x", "gyro_y", "gyro_z",
        "mag_x", "mag_y", "mag_z"
    ])
    csv_file.flush()


def cleanup_resources() -> None:
    global csv_file

    try:
        GPIO.cleanup()
    except Exception:
        pass

    try:
        if csv_file is not None:
            csv_file.flush()
            csv_file.close()
    except Exception:
        pass


# ============================================================
# Messlogik
# ============================================================

def read_one_measurement() -> MeasurementPoint:
    accel_x, accel_y, accel_z = sensor.acceleration
    gyro_x, gyro_y, gyro_z = sensor.gyro
    mag_x, mag_y, mag_z = sensor.magnetic
    timestamp = time.time()

    point = MeasurementPoint(
        timestamp_unix=timestamp,
        timestamp_iso=unix_to_iso(timestamp),
        accel_x=float(accel_x),
        accel_y=float(accel_y),
        accel_z=float(accel_z),
        gyro_x=float(gyro_x),
        gyro_y=float(gyro_y),
        gyro_z=float(gyro_z),
        mag_x=float(mag_x),
        mag_y=float(mag_y),
        mag_z=float(mag_z),
    )
    return point


def write_point_to_csv(point: MeasurementPoint) -> None:
    csv_writer.writerow([
        point.timestamp_unix,
        point.accel_x, point.accel_y, point.accel_z,
        point.gyro_x, point.gyro_y, point.gyro_z,
        point.mag_x, point.mag_y, point.mag_z
    ])
    csv_file.flush()


def measurement_loop() -> None:
    try:
        STATE.set_idle("Initialisiert")
        initialize_hardware()
        initialize_csv()
        STATE.clear_error()
    except Exception as exc:
        STATE.set_error(format_exception_message(exc))
        return

    try:
        while True:
            try:
                pin_state = GPIO.input(SWITCH_GPIO)

                if pin_state == GPIO.LOW:
                    point = read_one_measurement()
                    write_point_to_csv(point)
                    STATE.add_point(point)
                    STATE.clear_error()
                else:
                    STATE.set_idle("Nicht am Messen")

                time.sleep(MEASUREMENT_INTERVAL_SECONDS)

            except OSError as exc:
                # Hier landen oft echte I/O-Probleme
                STATE.set_error(format_exception_message(exc))
                time.sleep(1.0)

            except Exception as exc:
                STATE.set_error(format_exception_message(exc))
                time.sleep(1.0)

    finally:
        cleanup_resources()


# ============================================================
# Webserver / API
# ============================================================

@APP.route("/api/state")
def api_state():
    return jsonify(STATE.snapshot())


@APP.route("/")
def index():
    html = f"""<!DOCTYPE html>
<html lang="de">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover" />
  <title>Live-Messdaten Raspberry Pi</title>
  <style>
    :root {{
      --bg: #0b1220;
      --panel: rgba(17, 24, 39, 0.92);
      --text: #e5e7eb;
      --muted: #94a3b8;
      --ok: #16a34a;
      --warn: #dc2626;
      --idle: #d97706;
      --line: #60a5fa;
      --border: rgba(255,255,255,0.08);
      --shadow: 0 12px 30px rgba(0,0,0,0.25);
      --radius: 18px;
    }}

    * {{
      box-sizing: border-box;
    }}

    html, body {{
      margin: 0;
      padding: 0;
      background: linear-gradient(180deg, #08101d 0%, #111827 100%);
      color: var(--text);
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    }}

    body {{
      min-height: 100vh;
    }}

    .container {{
      width: min(1200px, calc(100vw - 24px));
      margin: 16px auto 32px auto;
    }}

    .card {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: var(--radius);
      box-shadow: var(--shadow);
      padding: 16px;
    }}

    .header {{
      position: sticky;
      top: 10px;
      z-index: 20;
      backdrop-filter: blur(10px);
    }}

    .header-row {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      gap: 12px;
      flex-wrap: wrap;
    }}

    h1, h2 {{
      margin: 0;
    }}

    h1 {{
      font-size: clamp(1.2rem, 2.5vw, 1.8rem);
    }}

    h2 {{
      font-size: 1.05rem;
      margin-bottom: 12px;
    }}

    .status-badge {{
      display: inline-flex;
      align-items: center;
      gap: 10px;
      border-radius: 999px;
      padding: 10px 14px;
      font-weight: 700;
      border: 1px solid var(--border);
      background: rgba(255,255,255,0.04);
    }}

    .status-dot {{
      width: 12px;
      height: 12px;
      border-radius: 50%;
      flex: 0 0 auto;
    }}

    .ok .status-dot {{
      background: var(--ok);
      box-shadow: 0 0 0 6px rgba(22,163,74,0.16);
    }}

    .warn .status-dot {{
      background: var(--warn);
      box-shadow: 0 0 0 6px rgba(220,38,38,0.16);
    }}

    .idle .status-dot {{
      background: var(--idle);
      box-shadow: 0 0 0 6px rgba(217,119,6,0.16);
    }}

    .meta {{
      margin-top: 12px;
      display: flex;
      flex-wrap: wrap;
      gap: 14px;
      color: var(--muted);
      font-size: 0.95rem;
    }}

    .error-box {{
      margin-top: 14px;
      display: none;
      padding: 14px;
      border-radius: 14px;
      border: 1px solid rgba(220,38,38,0.45);
      background: rgba(220,38,38,0.12);
      color: #fecaca;
      white-space: pre-wrap;
      word-break: break-word;
    }}

    .grid {{
      display: grid;
      grid-template-columns: 1.25fr 1fr;
      gap: 16px;
      margin-top: 16px;
    }}

    .canvas-wrap {{
      border: 1px solid var(--border);
      border-radius: 14px;
      padding: 8px;
      background: rgba(255,255,255,0.03);
    }}

    #plotCanvas {{
      width: 100%;
      height: 330px;
      display: block;
    }}

    .small {{
      color: var(--muted);
      font-size: 0.92rem;
      margin-top: 10px;
    }}

    .table-wrap {{
      overflow-x: auto;
    }}

    table {{
      width: 100%;
      border-collapse: collapse;
      min-width: 760px;
    }}

    th, td {{
      text-align: left;
      padding: 10px 8px;
      border-bottom: 1px solid var(--border);
      font-size: 0.93rem;
      font-variant-numeric: tabular-nums;
    }}

    th {{
      color: var(--muted);
      font-weight: 700;
    }}

    @media (max-width: 950px) {{
      .grid {{
        grid-template-columns: 1fr;
      }}
    }}

    @media (max-width: 640px) {{
      .container {{
        width: min(100vw - 12px, 1200px);
        margin: 10px auto 20px auto;
      }}

      .card {{
        padding: 14px;
        border-radius: 14px;
      }}

      #plotCanvas {{
        height: 260px;
      }}
    }}
  </style>
</head>
<body>
  <div class="container">
    <div class="card header">
      <div class="header-row">
        <h1>Live-Messdaten vom Raspberry Pi</h1>
        <div id="statusBadge" class="status-badge idle">
          <span class="status-dot"></span>
          <span id="statusText">Lade Status…</span>
        </div>
      </div>

      <div class="meta">
        <div>Aktualisierung: alle {API_REFRESH_INTERVAL_MS} ms</div>
        <div id="lastUpdateText">Letzter Messpunkt: —</div>
        <div id="pointCountText">Messpunkte im Speicher: —</div>
      </div>

      <div id="errorBox" class="error-box"></div>
    </div>

    <div class="grid">
      <div class="card">
        <h2>Liveplot: accel_z über echter Zeit</h2>
        <div class="canvas-wrap">
          <canvas id="plotCanvas"></canvas>
        </div>
        <div class="small" id="plotInfo">
          Die x-Achse verwendet die echten Zeitabstände zwischen den Messpunkten.
        </div>
      </div>

      <div class="card">
        <h2>Letzte {TABLE_POINTS} Messpunkte</h2>
        <div class="table-wrap">
          <table>
            <thead>
              <tr>
                <th>Zeit</th>
                <th>ax</th>
                <th>ay</th>
                <th>az</th>
                <th>gx</th>
                <th>gy</th>
                <th>gz</th>
                <th>mx</th>
                <th>my</th>
                <th>mz</th>
              </tr>
            </thead>
            <tbody id="tableBody">
              <tr><td colspan="10">Noch keine Daten</td></tr>
            </tbody>
          </table>
        </div>
      </div>
    </div>
  </div>

  <script>
    const REFRESH_MS = {API_REFRESH_INTERVAL_MS};
    const PLOT_WINDOW_SECONDS = {PLOT_WINDOW_SECONDS};

    const statusBadge = document.getElementById("statusBadge");
    const statusText = document.getElementById("statusText");
    const lastUpdateText = document.getElementById("lastUpdateText");
    const pointCountText = document.getElementById("pointCountText");
    const errorBox = document.getElementById("errorBox");
    const tableBody = document.getElementById("tableBody");
    const plotInfo = document.getElementById("plotInfo");

    const canvas = document.getElementById("plotCanvas");
    const ctx = canvas.getContext("2d");

    function formatNumber(x) {{
      if (x === null || x === undefined || Number.isNaN(x)) {{
        return "—";
      }}
      return Number(x).toFixed(4);
    }}

    function setStatus(isMeasuring, status, errorMessage) {{
      statusText.textContent = status || "Unbekannt";
      statusBadge.classList.remove("ok", "warn", "idle");

      if (errorMessage) {{
        statusBadge.classList.add("warn");
        errorBox.style.display = "block";
        errorBox.textContent = "Warnung / Fehler: " + errorMessage;
      }} else {{
        errorBox.style.display = "none";
        errorBox.textContent = "";
        if (isMeasuring) {{
          statusBadge.classList.add("ok");
        }} else {{
          statusBadge.classList.add("idle");
        }}
      }}
    }}

    function updateTable(points) {{
      if (!points || points.length === 0) {{
        tableBody.innerHTML = '<tr><td colspan="10">Noch keine Daten</td></tr>';
        return;
      }}

      const rows = points.slice().reverse().map(p => {{
        const t = String(p.timestamp_iso || "—").replace("T", " ");
        return `
          <tr>
            <td>${{t}}</td>
            <td>${{formatNumber(p.accel_x)}}</td>
            <td>${{formatNumber(p.accel_y)}}</td>
            <td>${{formatNumber(p.accel_z)}}</td>
            <td>${{formatNumber(p.gyro_x)}}</td>
            <td>${{formatNumber(p.gyro_y)}}</td>
            <td>${{formatNumber(p.gyro_z)}}</td>
            <td>${{formatNumber(p.mag_x)}}</td>
            <td>${{formatNumber(p.mag_y)}}</td>
            <td>${{formatNumber(p.mag_z)}}</td>
          </tr>
        `;
      }});

      tableBody.innerHTML = rows.join("");
    }}

    function resizeCanvas() {{
      const dpr = window.devicePixelRatio || 1;
      const rect = canvas.getBoundingClientRect();
      canvas.width = Math.max(300, Math.floor(rect.width * dpr));
      canvas.height = Math.max(220, Math.floor(rect.height * dpr));
      ctx.setTransform(1, 0, 0, 1, 0, 0);
      ctx.scale(dpr, dpr);
    }}

    function drawPlot(plotPoints) {{
      const rect = canvas.getBoundingClientRect();
      const width = rect.width;
      const height = rect.height;

      ctx.clearRect(0, 0, width, height);

      const pad = {{ left: 52, right: 18, top: 16, bottom: 36 }};
      const innerW = Math.max(10, width - pad.left - pad.right);
      const innerH = Math.max(10, height - pad.top - pad.bottom);

      ctx.fillStyle = "rgba(255,255,255,0.02)";
      ctx.fillRect(pad.left, pad.top, innerW, innerH);

      if (!plotPoints || plotPoints.length < 2) {{
        ctx.fillStyle = "#94a3b8";
        ctx.font = "14px system-ui, sans-serif";
        ctx.fillText("Noch nicht genug Daten für den Plot", pad.left + 10, pad.top + 22);
        return;
      }}

      const newestT = plotPoints[plotPoints.length - 1].timestamp_unix;
      const minAllowed = newestT - PLOT_WINDOW_SECONDS;
      const filtered = plotPoints.filter(p => p.timestamp_unix >= minAllowed);

      if (filtered.length < 2) {{
        ctx.fillStyle = "#94a3b8";
        ctx.font = "14px system-ui, sans-serif";
        ctx.fillText("Noch nicht genug Daten im gewählten Zeitfenster", pad.left + 10, pad.top + 22);
        return;
      }}

      const t0 = filtered[0].timestamp_unix;
      const xs = filtered.map(p => p.timestamp_unix - t0);
      const ys = filtered.map(p => Number(p.accel_z));

      let minX = Math.min(...xs);
      let maxX = Math.max(...xs);
      let minY = Math.min(...ys);
      let maxY = Math.max(...ys);

      if (maxX <= minX) {{
        maxX = minX + 1e-6;
      }}
      if (maxY <= minY) {{
        const eps = Math.abs(minY) * 0.05 + 1e-6;
        minY -= eps;
        maxY += eps;
      }}

      const yPad = Math.max((maxY - minY) * 0.10, 1e-6);
      minY -= yPad;
      maxY += yPad;

      ctx.strokeStyle = "rgba(148,163,184,0.18)";
      ctx.lineWidth = 1;

      const gridY = 4;
      for (let i = 0; i <= gridY; i++) {{
        const y = pad.top + (i / gridY) * innerH;
        ctx.beginPath();
        ctx.moveTo(pad.left, y);
        ctx.lineTo(pad.left + innerW, y);
        ctx.stroke();
      }}

      const gridX = 5;
      for (let i = 0; i <= gridX; i++) {{
        const x = pad.left + (i / gridX) * innerW;
        ctx.beginPath();
        ctx.moveTo(x, pad.top);
        ctx.lineTo(x, pad.top + innerH);
        ctx.stroke();
      }}

      ctx.strokeStyle = "rgba(229,231,235,0.45)";
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(pad.left, pad.top);
      ctx.lineTo(pad.left, pad.top + innerH);
      ctx.lineTo(pad.left + innerW, pad.top + innerH);
      ctx.stroke();

      ctx.fillStyle = "#94a3b8";
      ctx.font = "12px system-ui, sans-serif";

      for (let i = 0; i <= gridY; i++) {{
        const frac = 1 - i / gridY;
        const val = minY + frac * (maxY - minY);
        const y = pad.top + i / gridY * innerH;
        ctx.fillText(val.toFixed(3), 6, y + 4);
      }}

      for (let i = 0; i <= gridX; i++) {{
        const frac = i / gridX;
        const val = minX + frac * (maxX - minX);
        const x = pad.left + frac * innerW;
        ctx.fillText(val.toFixed(2) + " s", x - 10, pad.top + innerH + 18);
      }}

      ctx.strokeStyle = "#60a5fa";
      ctx.lineWidth = 2;
      ctx.beginPath();

      filtered.forEach((p, index) => {{
        const xVal = p.timestamp_unix - t0;
        const yVal = Number(p.accel_z);

        const x = pad.left + ((xVal - minX) / (maxX - minX)) * innerW;
        const y = pad.top + (1 - (yVal - minY) / (maxY - minY)) * innerH;

        if (index === 0) {{
          ctx.moveTo(x, y);
        }} else {{
          ctx.lineTo(x, y);
        }}
      }});

      ctx.stroke();

      const last = filtered[filtered.length - 1];
      const lx = pad.left + (((last.timestamp_unix - t0) - minX) / (maxX - minX)) * innerW;
      const ly = pad.top + (1 - ((Number(last.accel_z) - minY) / (maxY - minY))) * innerH;

      ctx.fillStyle = "#f8fafc";
      ctx.beginPath();
      ctx.arc(lx, ly, 3.5, 0, Math.PI * 2);
      ctx.fill();

      plotInfo.textContent =
        "Zeitfenster: " + (maxX - minX).toFixed(3) +
        " s | accel_z-Bereich: " + minY.toFixed(4) + " bis " + maxY.toFixed(4);
    }}

    async function fetchState() {{
      try {{
        const response = await fetch("/api/state", {{ cache: "no-store" }});
        if (!response.ok) {{
          throw new Error("HTTP " + response.status);
        }}

        const data = await response.json();

        setStatus(data.is_measuring, data.status_text, data.last_error);

        if (data.last_data_timestamp_unix && data.table_points.length > 0) {{
          const newest = data.table_points[data.table_points.length - 1].timestamp_iso;
          lastUpdateText.textContent = "Letzter Messpunkt: " + String(newest).replace("T", " ");
        }} else {{
          lastUpdateText.textContent = "Letzter Messpunkt: —";
        }}

        pointCountText.textContent = "Messpunkte im Speicher: " + (data.point_count ?? 0);

        updateTable(data.table_points || []);
        drawPlot(data.plot_points || []);
      }} catch (err) {{
        setStatus(false, "Verbindungsproblem", String(err));
      }}
    }}

    window.addEventListener("resize", () => {{
      resizeCanvas();
      fetchState();
    }});

    resizeCanvas();
    fetchState();
    setInterval(fetchState, REFRESH_MS);
  </script>
</body>
</html>
"""
    return Response(html, mimetype="text/html")


# ============================================================
# Start
# ============================================================

def main() -> None:
    thread = threading.Thread(target=measurement_loop, daemon=True)
    thread.start()
    APP.run(host=HOST, port=PORT, debug=False, threaded=True)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        cleanup_resources()