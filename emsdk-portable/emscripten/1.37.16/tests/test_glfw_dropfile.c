#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif
#define GLFW_INCLUDE_ES2
#include <GLFW/glfw3.h>

int result = 1;

GLFWwindow* g_window;

void render();

void render() {
  glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);
}

void on_file_drop(GLFWwindow *window, int count, const char **paths) {
  for (int i = 0; i < count; ++i) {
    printf("dropped file %s\n", paths[i]);

    FILE *fp = fopen(paths[i], "rb");
    if (!fp) {
        printf("failed to open %s\n", paths[i]);
        perror("fopen");
        result = 0;
        continue;
    }
    int c;
    long size = 0;
    bool dump = strstr(paths[i], ".txt") != 0;
    if (dump) printf("text file contents (first 100 bytes): ");
    while ((c = fgetc(fp)) != -1) {
        ++size;
        if (dump && size <= 100) putchar(c);
    }
    if (dump) putchar('\n');
    printf("read %ld bytes from %s\n", size, paths[i]);

    fclose(fp);

#ifdef __EMSCRIPTEN__
    // Emscripten copies the contents of the dropped file into the
    // in-browser filesystem. Delete after usage to free up memory.
    printf("unlinking %s\n", paths[i]);
    unlink(paths[i]);
#endif

  }
#ifdef REPORT_RESULT
  REPORT_RESULT();
#endif
}

int main() {
  if (!glfwInit())
  {
    result = 0;
    printf("Could not create window. Test failed.\n");      
#ifdef REPORT_RESULT
    REPORT_RESULT();
#endif      
    return -1;
  }
  glfwWindowHint(GLFW_RESIZABLE , 1);
  g_window = glfwCreateWindow(600, 450, "GLFW drop file", NULL, NULL);
  if (!g_window)
  {
    result = 0;
    printf("Could not create window. Test failed.\n");      
#ifdef REPORT_RESULT
    REPORT_RESULT();
#endif           
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(g_window);

  // Install callbacks
  glfwSetDropCallback(g_window, on_file_drop);

  // Main loop
  printf("Drag and drop a file from your desktop onto the green canvas.\n");
#ifdef __EMSCRIPTEN__
  emscripten_set_main_loop(render, 0, 1);
#else
  while (!glfwWindowShouldClose(g_window)) {
    render();
    glfwSwapBuffers(g_window);
    glfwPollEvents();
  }
#endif

  glfwTerminate();

  return 0;
}
