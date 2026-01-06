import logging
import torch
from dataclasses import dataclass
from networks.base_model import BaseModel
from normalizer import Normalizer


@dataclass
class TestResult:
    loss: float = 0
    accuracy: float = 0


class Tester():
    def __init__(self,
                 criterion: torch.nn.modules.loss._Loss,
                 testing_dataloader: torch.utils.data.DataLoader,
                 normalizer: Normalizer,
                 ):
        self.criterion = criterion
        self.testing_dataloader = testing_dataloader
        self.normalizer = normalizer
        self.logger = logging.getLogger(__name__)

    def test(self, network: BaseModel) -> TestResult:
        self.logger.info(f"Started testing network {network.get_name()}")

        # Load last parameters of the network
        network.load_model()
        network.eval()

        running_loss = 0
        correct = 0
        total = 0
        with torch.no_grad():
            for input, target in self.testing_dataloader:
                input = input.float()
                target = target.float()

                output = network(input)
                loss = self.criterion(output, target)
                running_loss += loss.item()

                _, target_denormalized = self.normalizer.denormalize(None, target)
                _, output_denormalized = self.normalizer.denormalize(None, output)

                predictions = output_denormalized.argmax(dim=1)
                correct += (predictions == target_denormalized).sum().item()
                total += target_denormalized.size(0)

        test_loss = running_loss / len(self.testing_dataloader)
        test_accuracy = correct / total

        result = TestResult(test_loss, test_accuracy)

        self.logger.info(f"Finished testing network {network.get_name()} with result\n{result}")

        return result


def test_all_models(testing_dataloader: torch.utils.data.DataLoader,
                    model_export_path: str,
                    normalizer: Normalizer) -> dict[str, TestResult]:
    # Loop over all defined networks and test them agains the testing data
    results = {}
    for network in BaseModel.subclasses:
        net = network(model_export_path)
        criterion = torch.nn.MSELoss()
        tester = Tester(criterion, testing_dataloader, normalizer)
        result = tester.test(net)
        results[net.get_name()] = result

    return results
