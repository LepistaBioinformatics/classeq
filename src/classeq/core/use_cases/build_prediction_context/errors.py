class UnexistingContext(Exception):
    def __init__(self, message: str) -> None:
        super().__init__(message)


class FailLoadingContext(Exception):
    def __init__(self, message: str) -> None:
        super().__init__(message)


class InvalidReferenceTreeFromContext(Exception):
    def __init__(self, message: str) -> None:
        super().__init__(message)


class UnexpectedErrorOnExecuteAdherenceTest(Exception):
    def __init__(self, message: str) -> None:
        super().__init__(message)
