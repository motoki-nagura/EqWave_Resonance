from dataclasses import dataclass
import numpy as np

@dataclass
class Params:
    day_in_sec: float
    N2: float
    n1_values: np.ndarray
    n2_values: np.ndarray
    n3_values: np.ndarray
    ω1: float
    ω2: float
    ω3: float
    k1_values: np.ndarray
    m1_guess: float
    k2_guesses: np.ndarray
    m2_guesses: np.ndarray
