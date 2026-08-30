class Book:
	def __init__(self, title, author, pages):
		self.title = title
		self.author = author
		self.pages = pages
	def summary(self):
		return f"{self.title}: the author is {self.author}, {self.pages}pages."
book = Book('Python', 'Jerry', 180)
print(book.summary())
