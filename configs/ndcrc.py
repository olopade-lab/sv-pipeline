from parsl.providers import CondorProvider

from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.utils import get_all_checkpoints

cores_per_slot = 8
worker_init = """

source ~/.bashrc
conda activate parsl
"""

scheduler_options = """
RequestMemory={}
""".format(5000 * cores_per_slot)

config = Config(
    executors=[
        HighThroughputExecutor(
            cores_per_worker=1,
            heartbeat_threshold=120,
            heartbeat_period=30,
            provider=CondorProvider(
                scheduler_options=scheduler_options,
                cores_per_slot=cores_per_slot,
                init_blocks=1,
                max_blocks=1,
                worker_init=worker_init,
            ),
        )
    ],
    checkpoint_mode='task_exit',
    checkpoint_files=get_all_checkpoints()
)
