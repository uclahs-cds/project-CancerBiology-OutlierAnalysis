import datetime
import tkinter as tk
from tkinter import ttk
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
import json


def format_time(seconds: int) -> str:
    """Format seconds into a time string."""
    h, m = divmod(seconds, 3600)
    m, s = divmod(m, 60)
    return f"{h:02}:{m:02}:{s:02}"


class CountdownTimerApp(tk.Tk):
    """A simple app to display timers."""
    def __init__(self):
        super().__init__()

        self.title("Multiple Countdown Timers")
        self.geometry("400x400")

        self.idle_timers = []  # List to hold timer objects
        self.active_timers = {}
        self.create_widgets()
        self.countdown_job = None

    def create_widgets(self):
        """Create timer display for up to 8 timers (without buttons)."""
        monospaced_font = ("Monaco", 12)
        for i in range(8):
            timer_frame = ttk.Frame(self)
            timer_frame.pack(pady=5)

            label_var = tk.StringVar()
            timer_var = tk.StringVar()
            # expected_var = tk.StringVar()

            ttk.Label(timer_frame, textvariable=label_var).pack(side="left", padx=5)

            timer_label = ttk.Label(
                timer_frame, textvariable=timer_var, font=monospaced_font
            )
            timer_label.pack(side="left", padx=5)

            timer = {
                "index": i,
                "label": label_var,
                "var": timer_var,
                "label_obj": timer_label,
                "expected": None,
            }
            self._make_available(timer)



    def start_timer(self, timer_name, seconds):
        """Set a timer programmatically using the given number of seconds."""
        self.idle_timers.sort(key=lambda x: x["index"], reverse=True)
        timer = self.idle_timers.pop()

        timer["label"].set(timer_name)
        # timer["expect_var"].set(format_time(seconds))

        timer["expected"] = datetime.datetime.now() + datetime.timedelta(
            seconds=seconds
        )

        self.active_timers[timer_name] = timer
        if not self.countdown_job:
            self.countdown_job = self.after(0, self.countdown)

    def stop_timer(self, timer_name):
        """Stop a timer by name."""
        self._make_available(self.active_timers.pop(timer_name))

    def _make_available(self, timer):
        """Stop the countdown timer."""
        timer["expected"] = None

        timer["label"].set("Idle")
        timer["var"].set(format_time(0))
        # timer["expect_var"].set(format_time(0))
        timer["label_obj"].config(foreground="black")
        self.idle_timers.append(timer)

    def countdown(self):
        """Countdown function for each timer."""
        for timer in self.active_timers.values():
            remaining = int(
                (timer["expected"] - datetime.datetime.now()).total_seconds()
            )
            timer["var"].set(format_time(abs(remaining)))

            if remaining < 0:
                timer["label_obj"].config(foreground="red")

        if self.active_timers:
            self.countdown_job = self.after(1000, self.countdown)
        else:
            self.countdown_job = None


# HTTP Server class to handle POST requests
class SimpleHTTPRequestHandler(BaseHTTPRequestHandler):
    def do_POST(self):
        content_length = int(self.headers["Content-Length"])
        post_data = self.rfile.read(content_length)

        try:
            data = json.loads(post_data)

            if "name" in data:
                if "seconds" in data:
                    timer_name = data["name"]
                    seconds = round(float(data["seconds"]))

                    # Set the timer via Tkinter's event loop to avoid threading issues
                    app.after(0, app.start_timer, timer_name, seconds)
                    self.send_response(200)
                    self.send_header("Content-type", "application/json")
                    self.end_headers()
                    response = {
                        "message": f"Timer {timer_name} started for {seconds} seconds"
                    }
                    self.wfile.write(json.dumps(response).encode())
                else:
                    timer_name = data["name"]

                    # Set the timer via Tkinter's event loop to avoid threading issues
                    app.after(0, app.stop_timer, timer_name)
                    self.send_response(200)
                    self.send_header("Content-type", "application/json")
                    self.end_headers()
                    response = {"message": f"Timer {timer_name} stopped"}
                    self.wfile.write(json.dumps(response).encode())
            else:
                self.send_error(400, "Invalid request body")
        except json.JSONDecodeError:
            self.send_error(400, "Invalid JSON")


# Function to run the HTTP server in a separate thread
def run_http_server():
    port = 65035
    server_address = ("", port)
    httpd = HTTPServer(server_address, SimpleHTTPRequestHandler)
    print(f"Starting HTTP server on port {port}")
    httpd.serve_forever()


if __name__ == "__main__":
    # Create and run the Tkinter app
    app = CountdownTimerApp()

    # Start the HTTP server in a separate thread
    server_thread = threading.Thread(target=run_http_server)
    server_thread.daemon = True  # Daemon thread will exit when the main program exits
    server_thread.start()

    # Start Tkinter's main loop
    app.mainloop()
