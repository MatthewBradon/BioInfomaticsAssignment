

str_input = str(input("Enter a string: "))

# Split the string into trigrams and if the last trigram is not complete just add it to the list
trigrams = [str_input[i:i+3] for i in range(0, len(str_input), 3)]
if len(str_input) % 3 != 0:
    trigrams[-1] = str_input[-(len(str_input) % 3):]


# Print out trigrams with a space in between
print(" ".join(trigrams))

