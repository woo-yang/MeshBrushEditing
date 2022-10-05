#pragma once
#include <GWindow.h>
#include "../common/mesh_type.h"

namespace gui {
	using namespace YRender;
	class Window :public GWindow {
	public:
		Window(std::string title, double width = 1200, double height = 600)
			:GWindow(title,width,height){};
	protected:
		void mouse_button_callback(GLFWwindow* window, int button, int action, int modifiers) override;
		void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos)override;
		void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)override;
		void keyboard_callback(GLFWwindow* window, int key, int scancode, int action, int mods)override;
	};

}