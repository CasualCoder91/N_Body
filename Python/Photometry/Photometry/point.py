from dataclasses import dataclass
import numpy as np

@dataclass
class Point:
    position: np.array([np.nan,np.nan])
    velocity: np.array([np.nan,np.nan])
    id: int #usually equal to id of the star
    magnitude: float

    def get_distance(self,pt2) -> float:
        return np.linalg.norm(self.position-pt2.position)

    def __array__(self) -> np.ndarray:
        return np.array([self.position[0], self.position[1], self.velocity[0],self.velocity[1]])


