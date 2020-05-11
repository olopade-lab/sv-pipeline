from parsl.providers import LocalProvider
from parsl.addresses import address_by_hostname
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.utils import get_all_checkpoints

config = Config(
    executors=[
        HighThroughputExecutor(
            address=address_by_hostname(),
            cores_per_worker=1,
            provider=LocalProvider(
                worker_init="conda activate parsl; module load singularity",
                init_blocks=1,
                max_blocks=1,
            ),
        )
    ],
    checkpoint_mode='task_exit',
    checkpoint_files=get_all_checkpoints()
)
