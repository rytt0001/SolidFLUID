extends Node

class StaticData:
	var totalParticles: int

# Called when the node enters the scene tree for the first time.
func _ready():
	var rd := RenderingServer.create_local_rendering_device()
	
	var shader_file := load("res://PhysicsComputer.glsl")
	var shader_spirv: RDShaderSPIRV = shader_file.get_spirv()
	var shader := rd.shader_create_from_spirv(shader_spirv)
	#Prepare our data. We use floats in the shader, so we need 32 bit.
	var input := PackedVector3Array([Vector3(1,0,0), Vector3(0,2,0), Vector3(0,0,3), Vector3(4,0,0), Vector3(0,5,0), Vector3(0,0,6), Vector3(7,0,0), Vector3(0,8,0), Vector3(0,0,9), Vector3(10,0,0)])
	var input_bytes := input.to_byte_array()
	var test := StaticData.new()
	test.totalParticles = 10;
	var staticInput := PackedInt32Array([test.totalParticles])
	var staticInput_bytes := staticInput.to_byte_array()
	#Create a storage buffer that can hold our float values.
	#Each float has 4 bytes (32 bit) so 10 x 4 = 40 bytes
	var staticData = rd.storage_buffer_create(staticInput_bytes.size(),staticInput_bytes)
	var particlesuniform = RDUniform.new()
	particlesuniform.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	particlesuniform.binding = 1
	
	particlesuniform.add_id(staticData)
	var buffer := rd.storage_buffer_create(input_bytes.size(), input_bytes)
	# Create a uniform to assign the buffer to the rendering device
	var uniform := RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	uniform.binding = 0 # this needs to match the "binding" in our shader file
	uniform.add_id(buffer)
	var uniform_set := rd.uniform_set_create([particlesuniform, uniform], shader, 0) # the last parameter (the 0) needs to match the "set" in our shader file
	# Create a compute pipeline
	var pipeline := rd.compute_pipeline_create(shader)
	var compute_list := rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, pipeline)
	rd.compute_list_bind_uniform_set(compute_list, uniform_set, 0)
	rd.compute_list_dispatch(compute_list, 5, 1, 1)
	rd.compute_list_end()
	rd.submit()
	rd.sync()
	# Read back the data from the buffer
	var output_bytes := rd.buffer_get_data(buffer)
	var output := output_bytes.to_float32_array()
	print("Input: ", input)
	print("Output: ", output)
	pass # Replace with function body.


# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta):
	pass
