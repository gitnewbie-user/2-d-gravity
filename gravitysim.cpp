// mass_sim_modern.cpp
// Build: C++17, GLFW, GLAD, OpenGL 3.3 Core
// 
// Windows (MSVC):
//   - Get GLAD (glad.h + glad.c for GL 3.3 core) and GLFW3
//   - cl /std:c++17 mass_sim_modern.cpp glad.c /Ipath\to\include /link glfw3.lib opengl32.lib gdi32.lib user32.lib shell32.lib
// Linux (g++):
//   g++ -std=c++17 mass_sim_modern.cpp glad.c -lglfw -ldl -lX11 -lpthread -o mass_sim
//
// This version replaces immediate-mode circles with a modern, efficient pipeline:
// - Instanced rendering of billboard quads
// - Fragment shader draws perfect circles (filled or outlined) with smooth edges
// - Vertex data sent once; per-object data via instance attributes
// - Adjustable smoothness (pixel-perfect, no triangle tesselation)

#include <cstdio>     // For printf, fprintf functions for console output
#include <cstdlib>    // For EXIT_FAILURE, EXIT_SUCCESS constants and exit functions
#include <cmath>      // For mathematical functions like sqrt, cbrt
#include <vector>     // For std::vector dynamic arrays
#include <string>     // For std::string class for text handling
#include <iostream>   // For std::cout, std::cerr stream operations

#include <glad/glad.h>  // OpenGL function loader - must come before GLFW
#include <GLFW/glfw3.h> // Cross-platform window and input handling library

using namespace std;    // Use standard namespace to avoid std:: prefix

// ------------------------------------------------------------
// Settings
// ------------------------------------------------------------
static int   gScreenW = 800;  // Window width in pixels
static int   gScreenH = 600;  // Window height in pixels
static const char* gTitle = "mass_sim (modern OpenGL)";  // Window title text

static const double kTargetDelta = 1.0 / 60.0;  // Target timestep: 60 Hz (16.67ms per frame)
static float gInitMass  = 50.0f;  // Default mass for newly created objects
static float gDensity   = 0.02f;  // Mass density in kg/m^3 (affects size calculation)

// Mouse interaction state variables
static bool   gMouseDown = false;  // True when left mouse button is pressed
static double gMouseX = 0.0, gMouseY = 0.0; // Mouse position in window pixels (origin top-left from GLFW)

// ------------------------------------------------------------
// Utility
// ------------------------------------------------------------
static float radiusFromMass(float mass, float density) {
    const float PI = 3.14159265359f;  // Mathematical constant pi
    float volume = mass / density; // Calculate volume from mass and density (V = m/ρ)
    // Calculate radius from volume using sphere formula: V = (4/3)πr³
    // Rearranged: r = cbrt(3V/(4π))
    return cbrt( (3.0f * volume) / (4.0f * PI) );
}

struct Object {
    // Position coordinates in screen pixels with origin at bottom-left (OpenGL convention)
    float x{400.0f};  // X position (horizontal, pixels from left edge)
    float y{300.0f};  // Y position (vertical, pixels from bottom edge)
    float vx{0.0f};   // X velocity component (pixels per second)
    float vy{0.0f};   // Y velocity component (pixels per second)
    float mass{gInitMass};  // Object mass (affects size and gravitational attraction)

    // State flags for object behavior
    bool initializing{false};  // True while user is creating/sizing the object
    bool launched{false};      // True after object has been launched with initial velocity

    // Visual properties - RGBA color components
    float r{1.0f}, g{0.2f}, b{0.2f}, a{1.0f};  // Red=1.0, Green=0.2, Blue=0.2, Alpha=1.0 (opaque red)

    // Rendering style options
    bool filled{true};              // True for filled circle, false for outline only
    float outlineThickness{2.0f};   // Thickness of outline in pixels (when not filled)

    // Calculate visual radius in pixels based on mass and density
    float radiusPx() const { return radiusFromMass(mass, gDensity); }

    // Apply acceleration to velocity (Newton's second law: a = F/m, so v += a*dt)
    void accelerate(float ax, float ay) { vx += ax; vy += ay; }
    
    // Update position using current velocity (Euler integration: x += v*dt)
    void update(float dt) { x += vx * dt; y += vy * dt; }

    // Handle collisions with window boundaries, applying energy loss
    void clampToBounds(float L, float R, float B, float T) {
        float rad = radiusPx();  // Get object's radius for collision detection
        
        // Left boundary collision
        if (x < L + rad) { x = L + rad; vx *= -0.8f; }  // Reflect X velocity with 80% energy retention
        
        // Right boundary collision  
        if (x > R - rad) { x = R - rad; vx *= -0.8f; }  // Reflect X velocity with 80% energy retention
        
        // Bottom boundary collision
        if (y < B + rad) { y = B + rad; vy *= -0.8f; }  // Reflect Y velocity with 80% energy retention
        
        // Top boundary collision
        if (y > T - rad) { y = T - rad; vy *= -0.8f; }  // Reflect Y velocity with 80% energy retention
    }
};

// ------------------------------------------------------------
// Shader helpers
// ------------------------------------------------------------
static GLuint compileShader(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);  // Create a new shader object of specified type
    glShaderSource(s, 1, &src, nullptr);  // Set the source code for the shader
    glCompileShader(s);  // Compile the shader source code
    
    // Check if compilation was successful
    GLint ok = 0; 
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);  // Get compilation status
    if (!ok) {  // If compilation failed
        // Get the length of the error log
        GLint len = 0; 
        glGetShaderiv(s, GL_INFO_LOG_LENGTH, &len);
        
        // Create a string to hold the error message
        string log(len, '\0');
        
        // Retrieve the actual error message
        glGetShaderInfoLog(s, len, nullptr, &log[0]);
        
        // Print error to stderr and exit program
        fprintf(stderr, "Shader compile error:\n%s\n", log.c_str());
        exit(EXIT_FAILURE);
    }
    return s;  // Return the compiled shader ID
}

static GLuint linkProgram(GLuint vs, GLuint fs) {
    GLuint p = glCreateProgram();  // Create a new shader program object
    glAttachShader(p, vs);         // Attach the vertex shader to the program
    glAttachShader(p, fs);         // Attach the fragment shader to the program
    glLinkProgram(p);              // Link the shaders together into a complete program
    
    // Check if linking was successful
    GLint ok = 0; 
    glGetProgramiv(p, GL_LINK_STATUS, &ok);  // Get the link status
    if (!ok) {  // If linking failed
        // Get the length of the error log
        GLint len = 0; 
        glGetProgramiv(p, GL_INFO_LOG_LENGTH, &len);
        
        // Create a string to hold the error message
        string log(len, '\0');
        
        // Retrieve the actual error message
        glGetProgramInfoLog(p, len, nullptr, &log[0]);
        
        // Print error to stderr and exit program
        fprintf(stderr, "Program link error:\n%s\n", log.c_str());
        exit(EXIT_FAILURE);
    }
    return p;  // Return the linked program ID
}

// ------------------------------------------------------------
// Circle billboard shaders (instanced)
// - Vertex: expands a unit quad around each instance center with correct pixel size
// - Fragment: computes distance to center; draws filled or outlined circle with AA
// ------------------------------------------------------------
static const char* kVS = R"GLSL(
#version 330 core
// Per-vertex attribute: corners of a unit square (shared by all instances)
layout (location = 0) in vec2 aQuadPos;      // Values: [-1,1]x[-1,1] unit square corners

// Per-instance attributes (different for each circle being drawn)
layout (location = 1) in vec2 iCenterPx;     // Circle center position in pixels (origin bottom-left)
layout (location = 2) in float iRadiusPx;    // Circle radius in pixels
layout (location = 3) in vec4  iColor;       // RGBA color values (0.0 to 1.0)
layout (location = 4) in float iFilled;      // >0.5 => filled circle, else outline only
layout (location = 5) in float iOutlinePx;   // Outline thickness in pixels (for non-filled circles)

// Uniform data (same for all vertices and instances in this draw call)
uniform vec2 uScreen; // Screen dimensions (width, height) in pixels

// Output data passed to fragment shader
out VS_OUT {
    vec2 localPos;     // Position within the quad [-1,1], used to compute distance from center
    vec4 color;        // Color to be used by fragment shader
    float radiusPx;    // Radius in pixels (passed through for fragment calculations)
    float filled;      // Fill flag (passed through for fragment shader)
    float outlinePx;   // Outline thickness (passed through for fragment shader)
} vs;

void main(){
    // Convert circle center from pixel coordinates to Normalized Device Coordinates (NDC)
    // Pixel coordinates: [0,screenWidth] x [0,screenHeight] (bottom-left origin)
    // NDC coordinates: [-1,1] x [-1,1] (OpenGL standard)
    vec2 centerNDC = vec2(
        (iCenterPx.x / uScreen.x) * 2.0 - 1.0,  // Convert X: pixels → [0,1] → [-1,1]
        (iCenterPx.y / uScreen.y) * 2.0 - 1.0   // Convert Y: pixels → [0,1] → [-1,1]
    );

    // Calculate how much to scale the unit quad in NDC space
    // We want the quad to be exactly 2*radius pixels in size
    vec2 scaleNDC = vec2( (iRadiusPx * 2.0) / uScreen.x,  // X scale factor in NDC
                          (iRadiusPx * 2.0) / uScreen.y ); // Y scale factor in NDC

    // Position this vertex by scaling the unit quad and centering it
    vec2 posNDC = centerNDC + aQuadPos * scaleNDC; // Final position in NDC

    // Set the final vertex position (OpenGL requires 4D homogeneous coordinates)
    gl_Position = vec4(posNDC, 0.0, 1.0);
    
    // Pass data to fragment shader
    vs.localPos   = aQuadPos;   // Local position within quad [-1,1] for distance calculations
    vs.color      = iColor;     // Circle color
    vs.radiusPx   = iRadiusPx;  // Circle radius in pixels
    vs.filled     = iFilled;    // Whether to fill or outline
    vs.outlinePx  = iOutlinePx; // Outline thickness
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core
// Input data from vertex shader (interpolated across the quad)
in VS_OUT {
    vec2 localPos;   // Position within quad [-1,1] x [-1,1] (interpolated from vertex shader)
    vec4 color;      // Circle color (RGBA)
    float radiusPx;  // Circle radius in pixels
    float filled;    // >0.5 means filled circle, ≤0.5 means outline only
    float outlinePx; // Thickness of outline in pixels (for non-filled circles)
} fs;

// Output: final pixel color
out vec4 FragColor;

void main(){
    // Calculate distance from fragment to circle center in pixels
    // localPos has length 1.0 at the circle edge, so multiply by radius to get pixel distance
    float distPx = length(fs.localPos) * fs.radiusPx;

    // Calculate anti-aliasing width (smooth transition zone around edges)
    // fwidth() gives the rate of change of distPx between adjacent pixels (~1 pixel)
    float edge = fs.radiusPx;  // The radius defines the circle edge
    float aa = fwidth(distPx); // Derivative-based anti-aliasing width

    float alpha = 0.0;  // Start with transparent pixel
    
    if (fs.filled > 0.5) {
        // FILLED CIRCLE: Inside circle = opaque (1.0), outside = transparent (0.0)
        // smoothstep creates smooth transition from 1 to 0 across the edge
        alpha = 1.0 - smoothstep(edge - aa, edge + aa, distPx);
    } else {
        // OUTLINE CIRCLE: Create a ring around the radius with specified thickness
        float halfT = max(fs.outlinePx * 0.5, 1.0);  // Half thickness, minimum 1 pixel
        float inner = edge - halfT;  // Inner edge of the ring
        float outer = edge + halfT;  // Outer edge of the ring
        
        // Calculate alpha for inside the ring (1.0 inside inner radius, 0.0 outside)
        float inside = 1.0 - smoothstep(inner - aa, inner + aa, distPx);
        
        // Calculate alpha for outside the ring (0.0 inside outer radius, 1.0 outside)  
        float outside = smoothstep(outer - aa, outer + aa, distPx);
        
        // Ring alpha = inside AND NOT outside
        alpha = clamp(inside - outside, 0.0, 1.0);
    }

    // Combine the base color with the calculated alpha
    FragColor = vec4(fs.color.rgb, fs.color.a * alpha);
    
    // Discard completely transparent pixels to improve performance
    if (FragColor.a <= 0.001) discard;
}
)GLSL";

// ------------------------------------------------------------
// GL resources
// ------------------------------------------------------------
static GLuint gProg = 0;         // Shader program ID (compiled vertex + fragment shaders)
static GLuint gVAO = 0;          // Vertex Array Object - stores vertex attribute configuration
static GLuint gVBOQuad = 0;      // Vertex Buffer Object for quad geometry (shared by all circles)
static GLuint gVBOInstances = 0; // Vertex Buffer Object for per-instance data (position, color, etc.)
static GLint  gLocScreen = -1;   // Uniform location for screen dimensions in shader

// Per-instance data layout - this struct matches the vertex shader instance attributes
struct Instance {
    float cx, cy;     // Circle center position in pixels
    float radius;     // Circle radius in pixels  
    float r, g, b, a; // RGBA color components (0.0 to 1.0)
    float filled;     // 1.0 for filled circle, 0.0 for outline only
    float outlinePx;  // Outline thickness in pixels
};

// ------------------------------------------------------------
// Globals
// ------------------------------------------------------------
static std::vector<Object> gObjs;  // Dynamic array storing all physics objects in the simulation

// ------------------------------------------------------------
// GLFW callbacks
// ------------------------------------------------------------
static void errorCallback(int code, const char* desc){
    // Called when GLFW encounters an error - print details to stderr
    fprintf(stderr, "GLFW error %d: %s\n", code, desc);
}

static void cursorPosCallback(GLFWwindow* wnd, double x, double y) {
    // Called whenever the mouse moves within the window
    gMouseX = x; gMouseY = y; // Store mouse position (y is top-left origin from GLFW)
    
    // Convert to OpenGL coordinate system (bottom-left origin) for display
    double mx = gMouseX;                  // X coordinate stays the same
    double my = gScreenH - gMouseY;       // Flip Y coordinate: top-left → bottom-left
    
    // Update window title to show current mouse coordinates
    string t = string(gTitle) + " | Mouse X: " + to_string((int)mx)
                  + " Y: " + to_string((int)my);
    glfwSetWindowTitle(wnd, t.c_str());   // Set the new window title
}

static void mouseButtonCallback(GLFWwindow* wnd, int button, int action, int mods) {
    // Called when mouse buttons are pressed or released
    if (button == GLFW_MOUSE_BUTTON_LEFT) {  // Only handle left mouse button
        if (action == GLFW_PRESS) {          // Mouse button pressed down
            gMouseDown = true;               // Set global mouse state flag
            
            // Convert mouse position from GLFW coordinates (top-left origin) to OpenGL (bottom-left)
            float x = (float)gMouseX;                    // X coordinate unchanged
            float y = (float)(gScreenH - gMouseY);       // Flip Y coordinate
            
            // Create a new object at the mouse position
            Object o; 
            o.x = x; o.y = y;           // Set position to mouse location
            o.vx = 0; o.vy = 0;         // Start with zero velocity
            o.mass = gInitMass;         // Use default initial mass
            o.initializing = true;      // Mark as being created (user can still modify)
            o.launched = false;         // Not yet launched with velocity
            
            gObjs.emplace_back(o);      // Add the new object to the simulation
            
        } else if (action == GLFW_RELEASE) { // Mouse button released
            gMouseDown = false;              // Clear global mouse state flag
            
            if (!gObjs.empty()) {            // If we have objects in the simulation
                Object& o = gObjs.back();    // Get reference to the most recently created object
                o.initializing = false;      // Object creation phase is complete
                o.launched = true;           // Mark as launched
                
                // Calculate launch velocity based on drag vector (mouse movement while held down)
                float tx = (float)gMouseX;                    // Current mouse X
                float ty = (float)(gScreenH - gMouseY);       // Current mouse Y (flipped)
                float dx = (tx - o.x);       // X displacement from object center to mouse
                float dy = (ty - o.y);       // Y displacement from object center to mouse
                float dist = sqrt(dx*dx + dy*dy);        // Total drag distance
                
                // Scale factor for velocity - longer drag = higher speed
                float k = 0.01f * dist;      // Tunable scaling factor
                
                // Set velocity opposite to drag direction (slingshot effect)
                o.vx = -dx * k;              // X velocity opposite to drag X
                o.vy = -dy * k;              // Y velocity opposite to drag Y
            }
        }
    }
}

// ------------------------------------------------------------
// Init GL
// ------------------------------------------------------------
static void initGL() {
    // Compile and link shaders
    GLuint vs = compileShader(GL_VERTEX_SHADER, kVS);    // Compile vertex shader from source
    GLuint fs = compileShader(GL_FRAGMENT_SHADER, kFS);  // Compile fragment shader from source
    gProg = linkProgram(vs, fs);                         // Link both shaders into a program
    glDeleteShader(vs); glDeleteShader(fs);              // Clean up individual shader objects (no longer needed)

    // Get the location of the uniform variable in the shader program
    gLocScreen = glGetUniformLocation(gProg, "uScreen");  // Find "uScreen" uniform for screen dimensions

    // Create geometry: a unit quad made of 2 triangles in local space [-1,1]
    // This will be instanced (reused) for every circle with different transformations
    const float quadVerts[12] = {
        -1.f,-1.f,  // Triangle 1: Bottom-left vertex
         1.f,-1.f,  // Triangle 1: Bottom-right vertex  
         1.f, 1.f,  // Triangle 1: Top-right vertex
        -1.f,-1.f,  // Triangle 2: Bottom-left vertex (shared)
         1.f, 1.f,  // Triangle 2: Top-right vertex (shared)
        -1.f, 1.f   // Triangle 2: Top-left vertex
    };

    // Create and configure Vertex Array Object (VAO)
    glGenVertexArrays(1, &gVAO);  // Generate one VAO
    glBindVertexArray(gVAO);      // Bind it as current (all subsequent vertex setup applies to this VAO)

    // Create and upload quad geometry to GPU
    glGenBuffers(1, &gVBOQuad);                                      // Generate buffer for quad vertices
    glBindBuffer(GL_ARRAY_BUFFER, gVBOQuad);                         // Bind as current array buffer
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);  // Upload vertex data (static = won't change)
    glEnableVertexAttribArray(0);                                    // Enable vertex attribute 0 (aQuadPos in shader)
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0);  // Configure: 2 floats per vertex, no offset

    // Create buffer for per-instance data (position, color, size, etc. for each circle)
    glGenBuffers(1, &gVBOInstances);                                           // Generate buffer for instance data
    glBindBuffer(GL_ARRAY_BUFFER, gVBOInstances);                              // Bind as current array buffer
    glBufferData(GL_ARRAY_BUFFER, 1024 * (GLsizeiptr)sizeof(Instance), nullptr, GL_DYNAMIC_DRAW);  // Allocate space for 1024 instances (dynamic = will change frequently)

    // Configure vertex attributes for per-instance data
    // Each Instance struct member becomes a vertex attribute sent to the shader
    size_t off = 0;  // Track byte offset within Instance struct
    
    // Attribute 1: Circle center position (2 floats: cx, cy)
    glEnableVertexAttribArray(1);                                              // Enable vertex attribute 1 (iCenterPx in shader)
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(Instance), (void*)off);  // 2 floats at current offset
    off += sizeof(float)*2;                                                   // Move offset past center (2 floats)

    // Attribute 2: Circle radius (1 float: radius)
    glEnableVertexAttribArray(2);                                              // Enable vertex attribute 2 (iRadiusPx in shader)
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(Instance), (void*)off);  // 1 float at current offset
    off += sizeof(float)*1;                                                   // Move offset past radius (1 float)

    // Attribute 3: Circle color (4 floats: r, g, b, a)
    glEnableVertexAttribArray(3);                                              // Enable vertex attribute 3 (iColor in shader)
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(Instance), (void*)off);  // 4 floats at current offset
    off += sizeof(float)*4;                                                   // Move offset past color (4 floats)

    // Attribute 4: Fill flag (1 float: filled)
    glEnableVertexAttribArray(4);                                              // Enable vertex attribute 4 (iFilled in shader)
    glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(Instance), (void*)off);  // 1 float at current offset
    off += sizeof(float)*1;                                                   // Move offset past filled flag (1 float)

    // Attribute 5: Outline thickness (1 float: outlinePx)
    glEnableVertexAttribArray(5);                                              // Enable vertex attribute 5 (iOutlinePx in shader)
    glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(Instance), (void*)off);  // 1 float at current offset

    // Set attribute divisors for instanced rendering
    // Divisor = 1 means the attribute advances once per instance (not per vertex)
    glVertexAttribDivisor(1, 1);  // Center position advances per instance
    glVertexAttribDivisor(2, 1);  // Radius advances per instance
    glVertexAttribDivisor(3, 1);  // Color advances per instance  
    glVertexAttribDivisor(4, 1);  // Fill flag advances per instance
    glVertexAttribDivisor(5, 1);  // Outline thickness advances per instance

    glBindVertexArray(0);  // Unbind VAO (clean state)

    // Configure OpenGL render state
    glEnable(GL_BLEND);                                 // Enable alpha blending for transparency
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // Standard alpha blending: src*alpha + dst*(1-alpha)
    glDisable(GL_DEPTH_TEST);                           // Disable depth testing (2D rendering)
    glClearColor(0,0,0,1);                              // Set clear color to black (R=0, G=0, B=0, A=1)
}

// ------------------------------------------------------------
// Render
// ------------------------------------------------------------
static void render() {
    glClear(GL_COLOR_BUFFER_BIT);  // Clear the screen to background color (black)

    // Convert physics objects to GPU-friendly instance data
    vector<Instance> instances;          // Create array to hold instance data
    instances.reserve(gObjs.size());          // Pre-allocate space for efficiency
    
    for (const Object& o : gObjs) {           // Loop through all physics objects
        Instance ins{};                       // Create new instance data structure
        ins.cx = o.x; ins.cy = o.y;          // Copy position from physics object
        ins.radius = o.radiusPx();            // Calculate visual radius from mass
        ins.r = o.r; ins.g = o.g; ins.b = o.b; ins.a = o.a;  // Copy RGBA color
        ins.filled = o.filled ? 1.0f : 0.0f;  // Convert boolean to float (1.0 = filled, 0.0 = outline)
        ins.outlinePx = o.outlineThickness;   // Copy outline thickness
        instances.push_back(ins);             // Add to instance array
    }

    // Set up GPU state for rendering
    glUseProgram(gProg);                                              // Activate our shader program
    glUniform2f(gLocScreen, (float)gScreenW, (float)gScreenH);       // Send screen dimensions to shader

    glBindVertexArray(gVAO);                                         // Bind our configured vertex array

    // Upload instance data to GPU
    glBindBuffer(GL_ARRAY_BUFFER, gVBOInstances);                    // Bind instance data buffer
    glBufferSubData(GL_ARRAY_BUFFER, 0,                              // Update buffer starting at offset 0
                    (GLsizeiptr)(instances.size()*sizeof(Instance)), // Size = number of instances * size per instance
                    instances.data());                               // Pointer to our instance data

    // Draw all circles in one efficient call
    // 6 vertices = 2 triangles per quad, instances.size() = number of circles to draw
    glDrawArraysInstanced(GL_TRIANGLES, 0, 6, (GLsizei)instances.size());

    // Clean up GPU state
    glBindVertexArray(0);  // Unbind vertex array
    glUseProgram(0);       // Unbind shader program
}

// ------------------------------------------------------------
// Physics (very simple — mimic original behavior)
// ------------------------------------------------------------
static void simulate(float dt) {
    // GRAVITATIONAL ATTRACTION: Each object attracts every other object
    // This is O(n²) complexity - every object affects every other object
    for (size_t i = 0; i < gObjs.size(); ++i) {        // For each object i
        for (size_t j = 0; j < gObjs.size(); ++j) {    // Check against every other object j
            if (i == j) continue;                       // Skip self-interaction
            
            // Calculate displacement vector from object i to object j
            float dx = gObjs[j].x - gObjs[i].x;         // X displacement
            float dy = gObjs[j].y - gObjs[i].y;         // Y displacement
            float d2 = dx*dx + dy*dy;                   // Distance squared (faster than sqrt)
            
            if (d2 > 1e-6f) {                           // Avoid division by zero for very close objects
                float d = sqrt(d2);                // Actual distance
                float ux = dx / d, uy = dy / d;         // Unit vector pointing from i to j
                
                // Apply simplified gravity force (not physically accurate, but stable for demo)
                // Real gravity would be: F = G*m1*m2/r², but we use constant force for simplicity
                gObjs[i].accelerate(ux * 20.0f * dt, uy * 20.0f * dt);  // Accelerate object i toward object j
            }
        }
    }

    // UPDATE OBJECT POSITIONS AND HANDLE SPECIAL STATES
    for (auto& o : gObjs) {                             // For each object in the simulation
        // OBJECT CREATION PHASE: While user is creating object (mouse held down)
        if (o.initializing) {
            // Allow user to increase mass by hovering mouse over the object
            float mx = (float)gMouseX;                               // Current mouse X (GLFW coordinates)
            float my = (float)(gScreenH - gMouseY);                  // Current mouse Y (converted to OpenGL coordinates)
            float rad = o.radiusPx();                                // Current object radius
            
            // Check if mouse is inside the object (simple box collision for performance)
            if (mx > o.x - rad && mx < o.x + rad && 
                my > o.y - rad && my < o.y + rad) {
                o.mass *= 1.01f;                                     // Slowly increase mass (1% per frame while hovering)
            }
            continue;                                                // Skip physics update while initializing
        }

        // PHYSICS INTEGRATION: Update position based on current velocity (Euler method)
        o.update(dt);                                                // Apply: position += velocity * deltaTime
        
        // BOUNDARY COLLISION: Keep objects within screen bounds with energy loss
        o.clampToBounds(0, (float)gScreenW, 0, (float)gScreenH);
    }

    // OBJECT-TO-OBJECT COLLISION DETECTION AND RESPONSE
    // Check every pair of objects for overlapping (O(n²) algorithm)
    for (size_t i = 0; i < gObjs.size(); ++i) {        // For each object i
        for (size_t j = i+1; j < gObjs.size(); ++j) {  // Check against objects j where j > i (avoid duplicate checks)
            // Calculate distance between object centers
            float dx = gObjs[j].x - gObjs[i].x;         // X displacement between centers
            float dy = gObjs[j].y - gObjs[i].y;         // Y displacement between centers
            float d = sqrt(dx*dx + dy*dy);         // Distance between centers
            float sumR = gObjs[i].radiusPx() + gObjs[j].radiusPx();  // Sum of both radii
            
            // Check for collision: distance < sum of radii means overlapping
            if (d < sumR && d > 1e-4f) {                // Collision detected (avoid division by zero)
                // COLLISION RESOLUTION: Separate overlapping objects
                float pen = sumR - d;                   // Penetration depth (how much they overlap)
                float ux = dx / d, uy = dy / d;         // Unit vector pointing from i to j
                
                // Push objects apart by half the penetration distance each
                gObjs[i].x -= ux * (pen*0.5f);         // Move object i away from j
                gObjs[i].y -= uy * (pen*0.5f);
                gObjs[j].x += ux * (pen*0.5f);         // Move object j away from i  
                gObjs[j].y += uy * (pen*0.5f);
                
                // ENERGY LOSS: Simple bounce with damping (not physically accurate but stable)
                gObjs[i].vy *= -0.8f;                   // Reverse and dampen Y velocity of object i
                gObjs[j].vy *= -0.8f;                   // Reverse and dampen Y velocity of object j
                // Note: Only Y velocity is affected - this is a simplified collision model
            }
        }
    }
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
int main() {
    // INITIALIZE GLFW LIBRARY
    glfwSetErrorCallback(errorCallback);            // Set callback function for GLFW errors
    if (!glfwInit()) {                              // Initialize GLFW library
        fprintf(stderr, "Failed to init GLFW\n"); 
        return EXIT_FAILURE; 
    }

    // CONFIGURE OPENGL CONTEXT
    // Request OpenGL 3.3 Core Profile (modern OpenGL, no deprecated features)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);  // OpenGL major version: 3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);  // OpenGL minor version: 3 (= OpenGL 3.3)
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // Core profile (no legacy OpenGL)

    // CREATE WINDOW
    GLFWwindow* wnd = glfwCreateWindow(gScreenW, gScreenH, gTitle, nullptr, nullptr);  // Create window with specified size and title
    if (!wnd) {                                     // Check if window creation failed
        fprintf(stderr, "Failed to create window\n"); 
        glfwTerminate();                            // Clean up GLFW
        return EXIT_FAILURE; 
    }
    glfwMakeContextCurrent(wnd);                    // Make this window's OpenGL context current
    glfwSwapInterval(1);                            // Enable V-Sync (limit to monitor refresh rate)

    // INITIALIZE OPENGL FUNCTION POINTERS
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {  // Load OpenGL functions using GLAD
        fprintf(stderr, "Failed to load GL via GLAD\n");
        return EXIT_FAILURE;
    }

    // REGISTER INPUT CALLBACKS
    glfwSetCursorPosCallback(wnd, cursorPosCallback);      // Register mouse movement callback
    glfwSetMouseButtonCallback(wnd, mouseButtonCallback);  // Register mouse button callback

    // INITIALIZE GRAPHICS SYSTEM
    initGL();                                       // Set up shaders, buffers, and OpenGL state

    // INITIALIZE SIMULATION
    gObjs.push_back(Object{});                      // Add one initial object to the simulation (like original)

    // MAIN LOOP TIMING VARIABLES
    double last = glfwGetTime();                    // Time of last frame
    double acc = 0.0;                               // Accumulator for fixed timestep

    // MAIN GAME LOOP
    while (!glfwWindowShouldClose(wnd)) {           // Continue until user closes window
        // CALCULATE FRAME TIME
        double now = glfwGetTime();                 // Current time
        double dt = now - last;                     // Delta time since last frame
        last = now;                                 // Update last frame time
        acc += dt;                                  // Accumulate time for physics

        // FIXED TIMESTEP PHYSICS UPDATE
        // Run physics at consistent rate regardless of rendering framerate
        while (acc >= kTargetDelta) {               // While we have enough accumulated time
            simulate((float)kTargetDelta);          // Run one physics step with fixed timestep
            acc -= kTargetDelta;                    // Subtract the time we just simulated
        }

        // RENDER FRAME
        render();                                   // Draw all objects to screen
        glfwSwapBuffers(wnd);                       // Swap front/back buffers (display the frame)
        glfwPollEvents();                           // Process window events (input, resize, etc.)
    }

    // CLEANUP
    glfwDestroyWindow(wnd);                         // Destroy the window
    glfwTerminate();                                // Clean up GLFW library
    return EXIT_SUCCESS;                            // Program completed successfully
}
