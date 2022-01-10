# Saves to the next 3-digit natural number available in the directory 
# this script is located within.

# import the necessary packages
import pyautogui
import os
import win32api
from PIL import ImageGrab
from functools import partial
ImageGrab.grab = partial(ImageGrab.grab, all_screens=True)

def couldbeint(v):
    try:
        int(v)
        return True
    except ValueError:
        return False

def maxA(arr):
    try:
        result = max(arr)
        return result if isinstance(result, int) else 0
    except ValueError:
        return 0

def get_next_num(path):
    # gets 3-wide
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    files = set(map(int, filter(couldbeint, map(lambda x: x[:3], files))))
    next_num = '{:>03}'.format(str(maxA(files)+1))
    return next_num

def get_region(A, B):
    # returns (left, top, width, height)
    X = [min(A[i], B[i]) for i in range(2)]
    Y = [max(A[i], B[i]) for i in range(2)]
    X.extend([Y[i]-X[i] for i in range(2)])
    return X

def main():
    print('Press Ctrl-C to quit.\n')
    print('First two left click is ignored to allow moving of terminal.')
    print('Screen specifications:')
    import ctypes
    from time import sleep
    user32 = ctypes.windll.user32
    screensize = user32.GetSystemMetrics(0), user32.GetSystemMetrics(1)
    print(f"{screensize = }")
    coord1 = [0, 0]
    state_left = win32api.GetKeyState(0x01)  # Left button down = 0 or 1. Button up = -127 or -128
    state_right = win32api.GetKeyState(0x02)  # Right button down = 0 or 1. Button up = -127 or -128
    
    try:
        flag = False
        i = -1
        while True:
            # prime mouse leys
            a = win32api.GetKeyState(0x01)
            b = win32api.GetKeyState(0x02)
            x, y = pyautogui.position()
            positionStr = 'X: ' + str(x).rjust(4) + ' Y: ' + str(y).rjust(4)
            print(positionStr, end='')
            print('\b' * len(positionStr), end='', flush=True)
            if a != state_left: # on mouse down
                state_left = a 
                if a < 0:
                    i += 1
                    if not flag:
                        # get coord 1
                        coord1 = pyautogui.position()
                        print(f"{coord1 = }")
                        
                    elif flag:
                        # get coord 2; and save image
                        coord2 = pyautogui.position()
                        print(f"{coord2 = }")
                        my_region = get_region(coord1, coord2)
                        print(f"xywh = {my_region}")
                        image_name = os.path.join(os.path.dirname(__file__),
                            get_next_num(os.path.dirname(__file__)) + '.png')
                        if i > 1:
                            pyautogui.screenshot(image_name, region= my_region)
                            print(f'Image has been saved to {image_name}')

                        # empty coord
                        coord1, coord2 = [0,0], [0,0]
                    flag = not flag 
                    
            elif b != state_right:
                print('Right-click, now terminating.')
                sleep(4)
                return
    except KeyboardInterrupt:
        print('Keyboard interrupt, now terminating.')
        sleep(4)
        return

# drive function
if __name__=="__main__":
    main()
    # get_next_num(r'D:\NUS\Semester 7\PC4199 Honours Project\20220114 Smarter Every Friday Presentation 2\images')  