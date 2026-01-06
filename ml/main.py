import os
import logging
from datetime import datetime
from data_loader import load_data
from trainer import train_all_models
from tester import test_all_models
from reporter import Reporter
from normalizer import Normalizer
from hyperparameters import Hyperparameters


def create_run_dir():
    # Create run folder
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_path = f"assets/runs/run-{timestamp}/"
    reports_export_path = os.path.join(run_path, "reports/")
    models_export_path = os.path.join(run_path, "models/")
    os.makedirs(run_path, exist_ok=True)
    os.makedirs(models_export_path, exist_ok=True)
    os.makedirs(reports_export_path, exist_ok=True)

    return (run_path, reports_export_path, models_export_path)


def main():
    range_info_path = "./data/min_max.yaml"
    (run_path, reports_export_path, models_export_path) = create_run_dir()

    # Configure logging
    logging.basicConfig(
        filename=os.path.join(run_path, "log.log"),
        level=logging.INFO,
        filemode="w")

    normalizer = Normalizer(range_info_path)
    reporter = Reporter(reports_export_path)
    (training_dataloader,
     validation_dataloader,
     testing_dataloader) = load_data(normalizer)

    train_results = train_all_models(training_dataloader,
                                     validation_dataloader,
                                     models_export_path)
    reporter.create_training_report(train_results)

    # Create reports
    test_results = test_all_models(testing_dataloader,
                                   models_export_path,
                                   normalizer)
    reporter.create_testing_report(test_results)


main()
