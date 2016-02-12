from Tkinter import *
import tkMessageBox
import threading
import time

# function to run a function in a new thread; intended to run clustal omega
def new_thread(tool_fn):
	xyz = threading.Thread(target=tool_fn, args= ())
	xyz.start()
	while xyz.isAlive():
		app.update()
	xyz.join()

app = Tk()
app.title("DemoLimo")
app.geometry("100x50")

def browse_fn():
	try:
		a= b+1
	except Exception as e:
		print 'Error that occured: %s' %e
	time.sleep(1)
	print 'asdf'
	tkMessageBox.showerror('Error', "Error occured!")


button_seq = Button(app, text = 'Run now', command=lambda: new_thread(browse_fn))
button_seq.pack()

app.mainloop()