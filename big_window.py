from interface import Interface

if __name__ == "__main__":
    interface = Interface(windowsize=(18, 12))
    interface.main_window.protocol("WM_DELETE_WINDOW", interface.on_closing)
    interface.main_window.mainloop()
