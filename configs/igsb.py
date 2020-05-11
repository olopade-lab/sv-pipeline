from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.addresses import address_by_hostname
from parsl.executors import HighThroughputExecutor
from parsl.monitoring.monitoring import MonitoringHub
from parsl.utils import get_all_checkpoints

config = Config(
    executors=[
        HighThroughputExecutor(
            worker_debug=True,
            max_workers=5,
            address=address_by_hostname(),
            provider=SlurmProvider(
                'daenerys',
                nodes_per_block=1,
                init_blocks=15,
                max_blocks=15,
                scheduler_options='#SBATCH --exclude=kg15-11', # docker: Got permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Post http://%2Fvar%2Frun%2Fdocker.sock/v1.37/containers/create: dial unix /var/run/docker.sock: connect: permission denied.
                # worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages',
                # worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages; docker pull olopadelab/polyfuse',
                worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages; docker load -i /cephfs/users/annawoodard/polyfuse/docker/polyfuse.tar',
                walltime='48:00:00'
            ),
        )
    ],
    checkpoint_mode='task_exit',
    checkpoint_files=get_all_checkpoints(),
   # monitoring=MonitoringHub(
   #     hub_address=address_by_hostname(),
   #     hub_port=55055,
   #     monitoring_debug=False,
   #     resource_monitoring_interval=10,
   # ),
   # retries=2,
)
