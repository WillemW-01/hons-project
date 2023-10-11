class Set:
    def __init__(self, items=None):
        self.bag = []
        if items is not None:
            self.add_items(items)

    def add(self, item):
        if item not in self.bag:
            self.bag.append(item)

    def add_items(self, items):
        for item in items:
            self.add(item)

    def index(self, value):
        try:
            return self.bag.index(value)
        except:
            return None

    def __iter__(self):
        return iter(self.bag)

    def __getitem__(self, index):
        return self.bag[index]

    def __setitem__(self, item, value):
        self.bag[item] = value

    def __len__(self):
        return len(self.bag)

    def __str__(self):
        return str(self.bag)
