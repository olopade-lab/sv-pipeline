from parsl.providers import SlurmProvider
from parsl.monitoring.monitoring import MonitoringHub
from parsl.addresses import address_by_hostname
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.utils import get_all_checkpoints

worker_init = """

module load singularity
conda activate parsl
"""

scheduler_options = """
#SBATCH --mem-per-cpu=40G
#SBATCH --cpus-per-task=1
"""

config = Config(
    executors=[
        HighThroughputExecutor(
            'htex',
            # mem_per_worker=30,
            # cores_per_worker=28,
            max_workers=1,
            heartbeat_threshold=360,
            heartbeat_period=30,
            worker_debug=True,
            address='10.50.220.71',
            provider=SlurmProvider(
                'broadwl',
                walltime='36:00:00',
                exclusive=False,
                scheduler_options=scheduler_options,
                worker_init=worker_init,
                init_blocks=70,
                max_blocks=70
            ),
        )
    ],
   monitoring=MonitoringHub(
       hub_address=address_by_hostname(),
       hub_port=55056,
       monitoring_debug=False,
       resource_monitoring_interval=60,
   ),
    checkpoint_mode='task_exit',
    checkpoint_files=get_all_checkpoints(),
    # retries=2
)
