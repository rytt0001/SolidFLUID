#[compute]
#version 450

// Invocations in the (x, y, z) dimension
layout(local_size_x = 2, local_size_y = 1, local_size_z = 1) in;


layout(set = 0, binding = 1, std430) restrict buffer StaticData {
    int TotalParticles;
} static_data;
// A binding to the buffer we create in our script
layout(set = 0, binding = 0, std430) restrict buffer MyDataBuffer {
    vec3 posData[];
}
my_data_buffer;


// The code we want to execute in each invocation
void main() {
    // gl_GlobalInvocationID.x uniquely identifies this invocation across all work groups
    my_data_buffer.posData[gl_GlobalInvocationID.x].x *= static_data.TotalParticles;
    my_data_buffer.posData[gl_GlobalInvocationID.x].y *= static_data.TotalParticles;
    my_data_buffer.posData[gl_GlobalInvocationID.x].z *= static_data.TotalParticles;
}