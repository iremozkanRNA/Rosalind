a = 324
b = 24
c = a - b
print('a - b is', c)

# Problem
# Given: Two positive integers a and b , each less than 1000.

# Return: The integer corresponding to the square of the hypotenuse of the right triangle whose legs have lengths a and b.

# Answer
x = 965
y = 886
hypotenuse_squared =  x**2 + y**2
print('The square of the hypotenuse is', hypotenuse_squared)

#String and Lists 
# Any item in a list can be accessed by its index, or the number that indicates its place in the list. For example, try running the following code:


tea_party = ['March Hare', 'Hatter', 'Dormouse', 'Alice']
print(tea_party[2])

tea_party[1] = 'Cheshire Cat'
print(tea_party)

tea_party.append('Jabberwocky')
print(tea_party)
# If you need to obtain only some of a list, you can use the notation list_name[a:b] to get only those from index a up to but not including index b. 
# For example, tea_party[1:3] returns Cheshire Cat, Dormouse, not Cheshire Cat, Dormouse, Alice. This process is called "list slicing".

# Problem
# Given: A string s of length at most 200 letters and four integers a, b, c and d.

# Return: The slice of this string from indices a through b and c through d (with space in between), inclusively. 
# In other words, we should include elements s[b]and s[d]in our slice.
# Problem
# Given: A string s of length at most 200 letters and four integers a, b, c and d.

# Return: The slice of this string from indices a through b and c through d (with space in between), inclusively. 
# In other words, we should include elements s[b]and s[d]in our slice.

train = "19EQYxl3lmh40Zq486z1STg4trRSIYJvCx1ff0rNTeImRnqbtBF8DB8P4VzL14Kinosternonu8ULTIK46vLs62xFGTqI32Gx55ft5ZTe4KGT5VTaSmWX1bOky68nu5TxBapi9Cc9cdFjIWnNlpalternaOiuOM."

print(train[62:73] +" " + train[147:154])
