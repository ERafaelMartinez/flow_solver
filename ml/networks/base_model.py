import os
import logging
import torch
import torch.nn as nn
from hyperparameters import Hyperparameters
from abc import ABC, abstractmethod


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
