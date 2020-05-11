from parsl.config import Config
from parsl.providers import SlurmProvider
from parsl.addresses import address_by_hostname
from parsl.executors import HighThroughputExecutor

config = Config(
    executors=[
        HighThroughputExecutor(
            worker_debug=True,
            max_workers=7,
            address=address_by_hostname(),
            provider=SlurmProvider(
                'daenerys',
                nodes_per_block=1,
                init_blocks=15,
                max_blocks=15,
                # worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages',
                # worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages; docker pull olopadelab/polyfuse',
                worker_init='docker stop $(docker ps -aq); export PYTHONPATH=$PYTHONPATH:/cephfs/users/annawoodard/.local/lib/python3.7/site-packages; docker load -i /cephfs/users/annawoodard/polyfuse/docker/polyfuse.tar',
                walltime='48:00:00'
            ),
        )
    ]
)
