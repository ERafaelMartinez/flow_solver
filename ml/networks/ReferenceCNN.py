from networks.base_model import BaseModel
from torch import Tensor
import torch.nn as nn
from hyperparameters import Hyperparameters


class ReferenceCNN(BaseModel):
    """
    Class which implements the model architecture of the
    flow simulator accelerator. It is based on nested
    2D convolutional layers.
    """

    @property
    def params(self) -> Hyperparameters:
        return Hyperparameters(epochs=100)

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
