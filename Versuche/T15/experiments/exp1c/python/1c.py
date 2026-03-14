import RPi.GPIO as GPIO
import time

GPIO.setmode(GPIO.BOARD)
GPIO.setwarnings(False)

led1 = 1
led2 = 18
led3 = 22

GPIO.setup(led1, GPIO.OUT)
GPIO.setup(led2, GPIO.OUT)
GPIO.setup(led3, GPIO.OUT)

def set_leds(a: bool, b: bool, c: bool):
    GPIO.output(led1, a)
    GPIO.output(led2, b)
    GPIO.output(led3, c)

with open("1c_output.txt", "w") as f:
    beispiel = [
        (0,0,0),
        (1,0,0),
        (0,1,0),
        (0,0,1),
        (1,1,0),
        (0,1,1),
        (1,0,1),
        (1,1,1)
    ]

    for i in range(5):
        for a,b,c in beispiel:
            set_leds(bool(a), bool(b), bool(c))

            timestamp = time.time()
            f.write(f"{timestamp}, {a}, {b}, {c}\n")
            time.sleep(1.5)   # längere Pause: besser sichtbar

            # kurze Aus-Phase zwischen Zuständen
            set_leds(False, False, False)
            time.sleep(0.3)
GPIO.cleanup()