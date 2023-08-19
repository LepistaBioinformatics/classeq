from multiprocessing import Process, Queue
from resource import RUSAGE_CHILDREN, getrusage
from typing import Any, Callable

from classeq.settings import LOGGER


def with_resource_monitoring(func: Callable[..., Any]) -> Callable[..., Any]:
    def queue_executor(
        queue: Queue,  # type: ignore
        func: Callable[..., Any],
        **kwargs: Any,
    ) -> None:
        queue.put(func(**kwargs))

    def wrapper(**kwargs: Any) -> None:
        try:
            LOGGER.debug(getrusage(RUSAGE_CHILDREN))

            queue: Any = Queue()

            proc = Process(
                target=queue_executor,
                args=(queue, func),
                kwargs=kwargs,
            )

            proc.start()
            proc.join()

            LOGGER.debug(getrusage(RUSAGE_CHILDREN))
            return queue.get()

            # return func(**kwargs)
        except Exception as exc:
            raise exc

    return wrapper
