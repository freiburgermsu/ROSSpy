a = input("WRITE SOME INPUT HERE: ")

print("I am printing your input: " + str(a))

f = open("test.txt", "w")
f.write("hello!")
f.close()