import json
import time
import warnings

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


def write_event_snapshot(df, save_snapshot_path: str, cols_to_keep, tree_name: str = "Events"):
    """
    Generated with Copilot GPT
    Write a small ROOT snapshot and a JSON file containing selected event information.

    Parameters
    - df: RDataFrame filtered to the events of interest
    - save_snapshot_path: base path (without extension) for output files
    - cols_to_keep: list of column names to write
    - tree_name: name of the ROOT tree to store (default: 'Events')

    The function will attempt to write both `{save_snapshot_path}.root` (via
    `df.Snapshot`) and `{save_snapshot_path}.json` (by calling `df.AsNumpy` and
    serializing entries). Warnings are issued on failure but the function
    attempts both operations independently.
    """

    # Attempt ROOT snapshot (best-effort)
    try:
        df.Snapshot(tree_name, f"{save_snapshot_path}.root", cols_to_keep)
    except Exception as e:
        warnings.warn(f"Failed to write ROOT snapshot ({save_snapshot_path}.root): {e}")

    # Attempt JSON export
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
                    # convert numpy scalar to python native
                    try:
                        if hasattr(v, 'item'):
                            v = v.item()
                    except Exception:
                        pass
                    # decode bytes if needed
                    if isinstance(v, (bytes, bytearray)):
                        try:
                            v = v.decode()
                        except Exception:
                            v = str(v)
                    entry[col] = v
                json_list.append(entry)

        with open(f"{save_snapshot_path}.json", "w") as jf:
            json.dump(json_list, jf, indent=2)
    except Exception as e:
        warnings.warn(f"Failed to write JSON snapshot ({save_snapshot_path}.json): {e}")