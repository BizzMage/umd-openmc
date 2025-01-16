import numpy as np

def in2cm(arr):
    """Converts list of measurements in inches to cm. Works on single numbers too.
    Paramters
    ---------
    arr : list of floats
        Measurements in inches
    Returns
    -------
    converted : list of floats
        Measurements converted to centimeters
    """
    try:
        return np.array([round(elem * 2.54, 8) for elem in arr])
    except:
        return round(arr*2.54, 8)

if __name__ == "__main__":
    print(in2cm((1,2, 3)))