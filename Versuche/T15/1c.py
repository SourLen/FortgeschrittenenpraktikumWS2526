import RPi.GPIO as GPIO
import time

GPIO.setmode(GPIO.BOARD)
GPIO.setwarnings(False)

led1 = 0000
led2 = 0000
led3 = 0000

GPIO.setup(led1,GPIO.OUT)
GPIO.setup(led2,GPIO.OUT)
GPIO.setup(led3,GPIO.OUT)

def set_leds(a: bool, b: bool, c: bool):

    GPIO.output(led1,a)
    GPIO.output(led2,b)
    GPIO.output(led3,c)

with open("led_data_1c.txt","w") as f:

    beispiel = [
        (0,0,0),
        (0,0,1),
        (0,1,0),
        (0,1,1),
        (1,0,0),
        (1,0,1),
        (1,1,0),
        (1,1,1)
    ]

    for i in range(5):

        for a,b,c in beispiel:
            set_leds(bool(a),bool(b),bool(c)) ### Das lowkey hässlich aber gut
            timestamp = time.time()
            f.write(f"{timestamp}: {a}, {b}, {c}\n")
            time.sleep(1)

GPIO.cleanup()