import pyautogui
import time

hours = ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']


for hour in hours:
    pyautogui.click(-1528,12)
    pyautogui.click(x=-972,y=550)
    pyautogui.press('backspace')
    pyautogui.press('backspace')
    pyautogui.typewrite(message = hour)

    pyautogui.click(x=-1156,y=704)
    time.sleep(3)

    pyautogui.click(button = 'right')
    time.sleep(1)
    pyautogui.click(x=-1067,y=465)
    time.sleep(2)
    pyautogui.typewrite(message = '20191215'+hour)
    pyautogui.click(x=-754,y=618)

    pyautogui.click(-1168,12)
    pyautogui.click(-1528,12)
    
    
