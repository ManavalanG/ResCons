from Tkinter import *
root = Tk()

global ad
ad = 'testing'

from multiprocessing import Process
def foo():
	global ad
	global top_var
	# ad = 'changed'
	# print 'yup'
	# print ad
	top_var = 'changed'
	print top_var
	return 79

def mult(tool):
	xyz = Process(target=tool, args= ())
	print '1'
	xyz.start()
	print '2'
	xyz.join()
	# xyz.terminate()


def top_win():
	top = Toplevel(root)
	top.geometry("")

	button_3 = Button(top, text='chnage global', command= lambda: mult(foo))
	button_3.pack()

	global top_var
	top_var = 34

button_1 = Button(root, text='chnage global', command= lambda: mult(foo))
button_1.pack()

button_2 = Button(root, text='open top', command= top_win)
button_2.pack()


if __name__ == "__main__":
	root.mainloop()