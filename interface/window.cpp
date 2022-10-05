#include "window.h"


namespace gui {

	void Window::mouse_button_callback(GLFWwindow* window, int button, int action, int modifiers) {
		std::cout << "test";
	}
	void Window::cursor_pos_callback(GLFWwindow* window, double xpos, double ypos) {

	}
	void Window::scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {

	}
	void Window::keyboard_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {

	}
}