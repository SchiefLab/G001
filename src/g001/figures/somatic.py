def format_xtick_labels(L):
    "Format the xtick labels to have the week and sample type"
    new_labels = []
    lookup = {
        "-4": "wk -4\nMBC",
        "3": "wk 3\nGC",
        "4": "wk 4\nMBC",
        "8": "wk 8\nMBC",
        "9": "wk 9\nPB",
        "10": "wk 10\nMBC",
        "11": "wk 11\nGC",
        "16": "wk 16\nMBC",
    }
    for x in L:
        text = x.get_text()
        new_text = lookup[text]
        new_labels.append(new_text)
    return new_labels
