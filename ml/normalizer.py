import os
import yaml
import torch


class Normalizer():
    def __init__(self, range_info_path: str):
        if not os.path.exists(range_info_path):
            raise FileNotFoundError("Range information file not found")

        # load the range information file
        with open(range_info_path, "r") as f:
            self.range_info = yaml.safe_load(f)

    def normalize(self, inputs: torch.Tensor, outputs: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        inputs_norm = inputs.clone() if inputs is not None else None
        outputs_norm = outputs.clone() if outputs is not None else None

        # normalize the data
        if inputs_norm is not None:
            inputs_norm = (
                (inputs - self.range_info["inputs"]["u"]["min"]) /
                (self.range_info["inputs"]["u"]["max"] -
                 self.range_info["inputs"]["u"]["min"])
            )
        if outputs_norm is not None:
            outputs_norm[:, 0, :, :] = (
                (outputs[:, 0, :, :] - self.range_info["labels"]["u"]["min"]) /
                (self.range_info["labels"]["u"]["max"] -
                 self.range_info["labels"]["u"]["min"])
            )
            outputs_norm[:, 1, :, :] = (
                (outputs[:, 1, :, :] - self.range_info["labels"]["v"]["min"]) /
                (self.range_info["labels"]["v"]["max"] -
                 self.range_info["labels"]["v"]["min"])
            )

        return inputs_norm, outputs_norm

    def denormalize(self, inputs: torch.Tensor, outputs: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        # clone the data to avoid modifying the original data
        inputs_denorm = inputs.clone() if inputs is not None else None
        outputs_denorm = outputs.clone() if outputs is not None else None

        # denormalize the data
        if inputs_denorm is not None:
            inputs_denorm = (
                inputs * (self.range_info["inputs"]["u"]["max"] - self.range_info["inputs"]["u"]["min"]) +
                self.range_info["inputs"]["u"]["min"]
            )
        if outputs_denorm is not None:
            outputs_denorm[:, 0, :, :] = (
                outputs[:, 0, :, :] * (self.range_info["labels"]["u"]["max"] - self.range_info["labels"]["u"]["min"]) +
                self.range_info["labels"]["u"]["min"]
            )
            outputs_denorm[:, 1, :, :] = (
                outputs[:, 1, :, :] * (self.range_info["labels"]["v"]["max"] - self.range_info["labels"]["v"]["min"]) +
                self.range_info["labels"]["v"]["min"]
            )

        return inputs_denorm, outputs_denorm
