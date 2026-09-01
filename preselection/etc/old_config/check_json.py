import json

# Open the file and load the data
with open('xsecs_13TeV_new.json', 'r') as file:
    data1 = json.load(file)


with open('xsecs_13TeV.json', 'r') as file:
    data2 = json.load(file)


common_keys = data1.keys() & data2.keys()
keys_only_in_2 = data2.keys() - data1.keys()
keys_only_in_1 = data1.keys() - data2.keys()

print("common_keys", sorted(common_keys))

print("only in the old+++++++++")

print("keys_only_in_2", sorted(keys_only_in_2))

print("only in the new+++++++++")

print("keys_only_in_1", sorted(keys_only_in_1))


print("============================")

# Open the file and load the data
with open('xsecs_13p6TeV_new.json', 'r') as file:
    data1 = json.load(file)


with open('xsecs_13p6TeV.json', 'r') as file:
    data2 = json.load(file)


common_keys = data1.keys() & data2.keys()
keys_only_in_2 = data2.keys() - data1.keys()
keys_only_in_1 = data1.keys() - data2.keys()

print("common_keys", sorted(common_keys))

print("only in the old+++++++++")

print("keys_only_in_2", sorted(keys_only_in_2))

print("only in the new+++++++++")

print("keys_only_in_1", sorted(keys_only_in_1))
