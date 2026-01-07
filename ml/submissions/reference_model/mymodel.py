from torch import Tensor
import torch.nn as nn
from dataclasses import dataclass
import os
import logging
import torch
from abc import ABC, abstractmethod


@dataclass
class Hyperparameters:
    learning_rate: float = 0.001
    epochs: int = 4000


class BaseModel(nn.Module, ABC):
    subclasses = []

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.subclasses.append(cls)

    def __init__(self, model_export_path):
        super().__init__()
        self.logger = logging.getLogger(__name__)
        self.model_export_path = model_export_path

    def get_name(self):
        return type(self).__name__

    @property
    @abstractmethod
    def params() -> Hyperparameters:
        pass

    def load_model(self) -> None:
        """Load model weights from binary file"""
        name = self.get_name()

        # check that model file exist in the given path
        if not os.path.exists(os.path.join(self.model_export_path, name + ".pt")):
            raise FileNotFoundError("Model file not found")

        self.logger.debug(f"Loading model {name} from {self.model_export_path}")

        self.load_state_dict(torch.load(os.path.join(self.model_export_path, name + ".pt")))

    def export_model(self) -> None:
        """Export model weights to binary file"""
        name = self.get_name()

        os.makedirs(self.model_export_path, exist_ok=True)

        self.logger.debug(f"Saving model {name} to {self.model_export_path}")

        torch.save(self.state_dict(), os.path.join(self.model_export_path, name + ".pt"))


class Model(BaseModel):
    """
    Class which implements the model architecture of the
    flow simulator accelerator. It is based on nested
    2D convolutional layers.
    """

    @property
    def params(self) -> Hyperparameters:
        return Hyperparameters()

    def __init__(self, model_export_path: str):
        super().__init__(model_export_path)

        # First convolutional layer:
        # Conv2d(in=1, out=16, kernel=7x7, padding="same", stride=1)
        self.conv1 = nn.Conv2d(in_channels=1, out_channels=16, kernel_size=7,
                               stride=1, padding='same')
        # Generic activation layer:
        self.relu = nn.ReLU()

        # Hidden layers:
        n_hidden = 5
        self.hidden = nn.ModuleList([
            nn.Conv2d(
                in_channels=16, out_channels=16, kernel_size=7,
                stride=1, padding='same'
            )
            for _ in range(n_hidden)
        ])

        # Output layer:
        self.output_layer = nn.Conv2d(
            in_channels=16, out_channels=2, kernel_size=7,
            stride=1, padding='same'
        )

    def forward(self, x: Tensor):
        """Evaluate the model on a given input"""
        # Apply first + hidden layers following of an
        # activation relu-layer
        x = self.conv1(x)
        x = self.relu(x)
        for layer in self.hidden:
            x = layer(x)
            x = self.relu(x)

        # Apply output layer and return
        return self.output_layer(x)


def init_my_model():
    model = Model("./")
    return model


init_my_model()
