# hyperspy-filtering

The application was created in connection with my engineer's thesis. The idea of the program was to make an easy and convinient tool for materials engineers that would help them analyse materials' microstructure images in an automated matter. The logic of the application is fully based on the code of Piotr Macioł, PhD. Eng.

## Installation

In order to make the application runnable, please follow below steps:

1. Download the latest release of Hyperspy-bundle from https://github.com/hyperspy/hyperspy-bundle/releases
2. Run the installer that will install the distribution in the current directory.
3. Open the newly created folder, open WinPython Command Prompt and run the commands below:
  1. `pip install seaborn`
  2. `pip install django`
  3. `pip install mpld3`
4. Move to the application's directory (precisely the directory with Manage.py file) and run `python manage.py runserver` command
