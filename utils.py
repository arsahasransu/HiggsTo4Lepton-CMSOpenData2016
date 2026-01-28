import json
import time
import warnings
import numpy as np

# Decorator to measure the execution time of a function
def time_eval(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"Execution time of function {func.__name__}: {execution_time:.6f} seconds")
        return result
    return wrapper


def write_event_snapshot(df, save_snapshot_path: str, cols_to_keep: list, tree_name: str = "Events"):

    def convert_to_serializable(obj):
        """Convert non-JSON-serializable objects to JSON-compatible types."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.integer, np.floating)):
            return obj.item()
        elif isinstance(obj, (bytes, bytearray)):
            try:
                return obj.decode()
            except Exception:
                return str(obj)
        elif type(obj).__name__ == 'RVec':
            try:
                return list(obj)
            except Exception:
                return str(obj)
        elif hasattr(obj, 'item'):
            try:
                return obj.item()
            except Exception:
                return str(obj)
        return obj

    # Attempt JSON export with RVec and ndarray support
    try:
        arrs = df.AsNumpy(cols_to_keep)
        if not arrs:
            json_list = []
        else:
            # All arrays should have same length
            length = len(next(iter(arrs.values())))
            json_list = []
            for i in range(length):
                entry = {}
                for col in cols_to_keep:
                    vals = arrs.get(col)
                    if vals is None:
                        entry[col] = None
                        continue
                    v = vals[i]
                    entry[col] = convert_to_serializable(v)
                json_list.append(entry)

        with open(f"{save_snapshot_path}.json", "w") as jf:
            json.dump(json_list, jf, indent=2)
        print(f"Successfully wrote event snapshot to {save_snapshot_path}.json")
    except Exception as e:
        warnings.warn(f"Failed to write JSON snapshot ({save_snapshot_path}.json): {e}")