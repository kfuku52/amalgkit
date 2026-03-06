class AmalgkitExit(Exception):
    def __init__(self, message='', exit_code=1, use_stderr=True):
        super().__init__(message)
        self.message = '' if message is None else str(message)
        self.exit_code = int(exit_code)
        self.use_stderr = bool(use_stderr)
