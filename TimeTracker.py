import datetime


# Class to simulate Java example for timing and timestamping
class TimeTracker:
    def __init__(self):
        # Start time
        self.start_time = datetime.datetime.now()

    def get_timestamp(self):
        # Current time formatted as 'yyyy-MM-dd-HHmmss'
        return self.start_time.strftime("%Y-%m-%d-%H%M%S")

    def calculate_duration(self):
        # End time
        end_time = datetime.datetime.now()
        # Duration in milliseconds
        duration = (end_time - self.start_time).total_seconds() * 1000

        # Convert milliseconds to minutes and seconds
        duration_minutes = int(duration // 60000)
        duration_seconds = int((duration % 60000) / 1000)

        return duration_minutes, duration_seconds

# Example usage (Commented out to prevent execution here):
# tracker = TimeTracker()
# timestamp = tracker.get_timestamp()
# minutes, seconds = tracker.calculate_duration()
# print(f"Timestamp: {timestamp}, Duration: {minutes} minutes and {seconds} seconds")
