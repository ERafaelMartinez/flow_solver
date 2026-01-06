from dataclasses import dataclass


@dataclass
class Hyperparameters:
    learning_rate: float = 0.001
    epochs: int = 4000
