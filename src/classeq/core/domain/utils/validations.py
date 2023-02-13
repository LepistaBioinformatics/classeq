def should_be_int(value: str) -> bool:
    try:
        int(value)
        return True
    except Exception:
        return False
