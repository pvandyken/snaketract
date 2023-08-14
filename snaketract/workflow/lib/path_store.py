import os
import hashlib
import more_itertools as itx
import json
from collections import UserDict
def get_hash(items: str):
    encoded = items.encode("utf-8")
    return hashlib.md5(encoded).hexdigest()


class PathStore(UserDict):
    def __init__(self):
        self.cache = ".snakemake/path-store.json"
        if os.path.exists(self.cache):
            with open(self.cache) as f:
                super().__init__(json.load(f))
        else:
            super().__init__()

    def __call__(self, wcard_name):
        def inner(wcards):
            try:
                path = self.data[wcards[wcard_name]]
            except KeyError as err:
                raise Exception(self.data) from err
            return path.format(**wcards)
        return inner
    
    def _to_hashmap(self, __mapping):
        keymap = {}
        hashmap = {}
        for key, value in __mapping.items():
            if isinstance(value, (dict, UserDict)):
                sub_keymap, sub_hashmap = self._to_hashmap(value)
                keymap[key] = sub_keymap
                hashmap.update(sub_hashmap)
                continue
            hashed = get_hash(str(value))
            keymap[key] = hashed
            hashmap[hashed] = value
        return keymap, hashmap

    def register_mapping(self, mapping):
        keymap, hashmap = self._to_hashmap(mapping)
        self.data.update(hashmap)
        self.save()
        return keymap

    def register(self, __item):
        hashed = get_hash(str(__item))
        # Make nested register calls safe
        if hashed in self.data:
            self.data.update({hashed: __item})
            return hashed
        self.data.update({hashed: __item})
        self.save()
        return hashed

    def save(self):
        data = {key: str(val) for key, val in self.data.items()}
        with open(self.cache, 'w') as f:
            json.dump(data, f)

    @classmethod
    def from_mapping(cls, mapping):
        inst = cls()
        return inst.register_mapping(mapping), inst


if __name__ == "__main__":
    paths, store = PathStore.from_mapping({
        "path one": {
            "sub prop": "Path/to/some/place",
            "another prop": "/my/special/path"
        },
        "path two": "Path/to/another/place"
    })

    print(store('ref')({'ref': paths['path two']}))
    print(store[paths['path one']['sub prop']])
    store.register_mapping({"an update": "/the updated path"})
    print(store)