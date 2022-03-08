#include "pch.h"
#include "Terrain.h"


Terrain::Terrain()
{
	m_terrainGeneratedToggle = false;
}


Terrain::~Terrain()
{
}

bool Terrain::Initialize(ID3D11Device* device, int terrainWidth, int terrainHeight)
{
	int index;
	float height = 0.0;
	bool result;

	// Save the dimensions of the terrain.
	m_terrainWidth = terrainWidth;
	m_terrainHeight = terrainHeight;

	m_frequency = m_terrainWidth / 20;
	m_amplitude = 3.0;
	m_wavelength = 1;

	// Create the structure to hold the terrain data.
	m_heightMap = new HeightMapType[m_terrainWidth * m_terrainHeight];
	if (!m_heightMap)
	{
		return false;
	}

	//this is how we calculate the texture coordinates first calculate the step size there will be between vertices. 
	float textureCoordinatesStep = 5.0f / m_terrainWidth;  //tile 5 times across the terrain. 
	// Initialise the data in the height map (flat).
	for (int j = 0; j<m_terrainHeight; j++)
	{
		for (int i = 0; i<m_terrainWidth; i++)
		{
			index = (m_terrainHeight * j) + i;

			m_heightMap[index].x = (float)i;
			m_heightMap[index].y = (float)height;
			m_heightMap[index].z = (float)j;

			//and use this step to calculate the texture coordinates for this point on the terrain.
			m_heightMap[index].u = (float)i * textureCoordinatesStep;
			m_heightMap[index].v = (float)j * textureCoordinatesStep;

		}
	}

	//even though we are generating a flat terrain, we still need to normalise it. 
	// Calculate the normals for the terrain data.
	result = CalculateNormals();
	if (!result)
	{
		return false;
	}

	// Initialize the vertex and index buffer that hold the geometry for the terrain.
	result = InitializeBuffers(device);
	if (!result)
	{
		return false;
	}

	
	return true;
}

void Terrain::Render(ID3D11DeviceContext * deviceContext)
{
	// Put the vertex and index buffers on the graphics pipeline to prepare them for drawing.
	RenderBuffers(deviceContext);
	deviceContext->DrawIndexed(m_indexCount, 0, 0);

	return;
}

bool Terrain::CalculateNormals()
{
	int i, j, index1, index2, index3, index, count;
	float vertex1[3], vertex2[3], vertex3[3], vector1[3], vector2[3], sum[3], length;
	DirectX::SimpleMath::Vector3* normals;
	

	// Create a temporary array to hold the un-normalized normal vectors.
	normals = new DirectX::SimpleMath::Vector3[(m_terrainHeight - 1) * (m_terrainWidth - 1)];
	if (!normals)
	{
		return false;
	}

	// Go through all the faces in the mesh and calculate their normals.
	for (j = 0; j<(m_terrainHeight - 1); j++)
	{
		for (i = 0; i<(m_terrainWidth - 1); i++)
		{
			index1 = (j * m_terrainHeight) + i;
			index2 = (j * m_terrainHeight) + (i + 1);
			index3 = ((j + 1) * m_terrainHeight) + i;

			// Get three vertices from the face.
			vertex1[0] = m_heightMap[index1].x;
			vertex1[1] = m_heightMap[index1].y;
			vertex1[2] = m_heightMap[index1].z;

			vertex2[0] = m_heightMap[index2].x;
			vertex2[1] = m_heightMap[index2].y;
			vertex2[2] = m_heightMap[index2].z;

			vertex3[0] = m_heightMap[index3].x;
			vertex3[1] = m_heightMap[index3].y;
			vertex3[2] = m_heightMap[index3].z;

			// Calculate the two vectors for this face.
			vector1[0] = vertex1[0] - vertex3[0];
			vector1[1] = vertex1[1] - vertex3[1];
			vector1[2] = vertex1[2] - vertex3[2];
			vector2[0] = vertex3[0] - vertex2[0];
			vector2[1] = vertex3[1] - vertex2[1];
			vector2[2] = vertex3[2] - vertex2[2];

			index = (j * (m_terrainHeight - 1)) + i;

			// Calculate the cross product of those two vectors to get the un-normalized value for this face normal.
			normals[index].x = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1]);
			normals[index].y = (vector1[2] * vector2[0]) - (vector1[0] * vector2[2]);
			normals[index].z = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0]);
		}
	}

	// Now go through all the vertices and take an average of each face normal 	
	// that the vertex touches to get the averaged normal for that vertex.
	for (j = 0; j<m_terrainHeight; j++)
	{
		for (i = 0; i<m_terrainWidth; i++)
		{
			// Initialize the sum.
			sum[0] = 0.0f;
			sum[1] = 0.0f;
			sum[2] = 0.0f;

			// Initialize the count.
			count = 0;

			// Bottom left face.
			if (((i - 1) >= 0) && ((j - 1) >= 0))
			{
				index = ((j - 1) * (m_terrainHeight - 1)) + (i - 1);

				sum[0] += normals[index].x;
				sum[1] += normals[index].y;
				sum[2] += normals[index].z;
				count++;
			}

			// Bottom right face.
			if ((i < (m_terrainWidth - 1)) && ((j - 1) >= 0))
			{
				index = ((j - 1) * (m_terrainHeight - 1)) + i;

				sum[0] += normals[index].x;
				sum[1] += normals[index].y;
				sum[2] += normals[index].z;
				count++;
			}

			// Upper left face.
			if (((i - 1) >= 0) && (j < (m_terrainHeight - 1)))
			{
				index = (j * (m_terrainHeight - 1)) + (i - 1);

				sum[0] += normals[index].x;
				sum[1] += normals[index].y;
				sum[2] += normals[index].z;
				count++;
			}

			// Upper right face.
			if ((i < (m_terrainWidth - 1)) && (j < (m_terrainHeight - 1)))
			{
				index = (j * (m_terrainHeight - 1)) + i;

				sum[0] += normals[index].x;
				sum[1] += normals[index].y;
				sum[2] += normals[index].z;
				count++;
			}

			// Take the average of the faces touching this vertex.
			sum[0] = (sum[0] / (float)count);
			sum[1] = (sum[1] / (float)count);
			sum[2] = (sum[2] / (float)count);

			// Calculate the length of this normal.
			length = sqrt((sum[0] * sum[0]) + (sum[1] * sum[1]) + (sum[2] * sum[2]));

			// Get an index to the vertex location in the height map array.
			index = (j * m_terrainHeight) + i;

			// Normalize the final shared normal for this vertex and store it in the height map array.
			m_heightMap[index].nx = (sum[0] / length);
			m_heightMap[index].ny = (sum[1] / length);
			m_heightMap[index].nz = (sum[2] / length);
		}
	}

	// Release the temporary normals.
	delete[] normals;
	normals = 0;

	return true;
}

void Terrain::Shutdown()
{
	// Release the index buffer.
	if (m_indexBuffer)
	{
		m_indexBuffer->Release();
		m_indexBuffer = 0;
	}

	// Release the vertex buffer.
	if (m_vertexBuffer)
	{
		m_vertexBuffer->Release();
		m_vertexBuffer = 0;
	}

	return;
}

bool Terrain::InitializeBuffers(ID3D11Device * device )
{
	VertexType* vertices;
	unsigned long* indices;
	D3D11_BUFFER_DESC vertexBufferDesc, indexBufferDesc;
	D3D11_SUBRESOURCE_DATA vertexData, indexData;
	HRESULT result;
	int index, i, j;
	int index1, index2, index3, index4; //geometric indices. 

	// Calculate the number of vertices in the terrain mesh.
	m_vertexCount = (m_terrainWidth - 1) * (m_terrainHeight - 1) * 6;

	// Set the index count to the same as the vertex count.
	m_indexCount = m_vertexCount;

	// Create the vertex array.
	vertices = new VertexType[m_vertexCount];
	if (!vertices)
	{
		return false;
	}

	// Create the index array.
	indices = new unsigned long[m_indexCount];
	if (!indices)
	{
		return false;
	}

	// Initialize the index to the vertex buffer.
	index = 0;

	for (j = 0; j<(m_terrainHeight - 1); j++)
	{
		for (i = 0; i<(m_terrainWidth - 1); i++)
		{
			index1 = (m_terrainHeight * j) + i;          // Bottom left.
			index2 = (m_terrainHeight * j) + (i + 1);      // Bottom right.
			index3 = (m_terrainHeight * (j + 1)) + i;      // Upper left.
			index4 = (m_terrainHeight * (j + 1)) + (i + 1);  // Upper right.

															 // Upper left.
			vertices[index].position = DirectX::SimpleMath::Vector3(m_heightMap[index3].x, m_heightMap[index3].y, m_heightMap[index3].z);
			vertices[index].normal = DirectX::SimpleMath::Vector3(m_heightMap[index3].nx, m_heightMap[index3].ny, m_heightMap[index3].nz);
			vertices[index].texture = DirectX::SimpleMath::Vector2(m_heightMap[index3].u, m_heightMap[index3].v);
			indices[index] = index;
			index++;

			// Upper right.
			vertices[index].position = DirectX::SimpleMath::Vector3(m_heightMap[index4].x, m_heightMap[index4].y, m_heightMap[index4].z);
			vertices[index].normal = DirectX::SimpleMath::Vector3(m_heightMap[index4].nx, m_heightMap[index4].ny, m_heightMap[index4].nz);
			vertices[index].texture = DirectX::SimpleMath::Vector2(m_heightMap[index4].u, m_heightMap[index4].v);
			indices[index] = index;
			index++;

			// Bottom left.
			vertices[index].position = DirectX::SimpleMath::Vector3(m_heightMap[index1].x, m_heightMap[index1].y, m_heightMap[index1].z);
			vertices[index].normal = DirectX::SimpleMath::Vector3(m_heightMap[index1].nx, m_heightMap[index1].ny, m_heightMap[index1].nz);
			vertices[index].texture = DirectX::SimpleMath::Vector2(m_heightMap[index1].u, m_heightMap[index1].v);
			indices[index] = index;
			index++;

			// Bottom left.
			vertices[index].position = DirectX::SimpleMath::Vector3(m_heightMap[index1].x, m_heightMap[index1].y, m_heightMap[index1].z);
			vertices[index].normal = DirectX::SimpleMath::Vector3(m_heightMap[index1].nx, m_heightMap[index1].ny, m_heightMap[index1].nz);
			vertices[index].texture = DirectX::SimpleMath::Vector2(m_heightMap[index1].u, m_heightMap[index1].v);
			indices[index] = index;
			index++;

			// Upper right.
			vertices[index].position = DirectX::SimpleMath::Vector3(m_heightMap[index4].x, m_heightMap[index4].y, m_heightMap[index4].z);
			vertices[index].normal = DirectX::SimpleMath::Vector3(m_heightMap[index4].nx, m_heightMap[index4].ny, m_heightMap[index4].nz);
			vertices[index].texture = DirectX::SimpleMath::Vector2(m_heightMap[index4].u, m_heightMap[index4].v);
			indices[index] = index;
			index++;

			// Bottom right.
			vertices[index].position = DirectX::SimpleMath::Vector3(m_heightMap[index2].x, m_heightMap[index2].y, m_heightMap[index2].z);
			vertices[index].normal = DirectX::SimpleMath::Vector3(m_heightMap[index2].nx, m_heightMap[index2].ny, m_heightMap[index2].nz);
			vertices[index].texture = DirectX::SimpleMath::Vector2(m_heightMap[index2].u, m_heightMap[index2].v);
			indices[index] = index;
			index++;
		}
	}

	// Set up the description of the static vertex buffer.
	vertexBufferDesc.Usage = D3D11_USAGE_DEFAULT;
	vertexBufferDesc.ByteWidth = sizeof(VertexType) * m_vertexCount;
	vertexBufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	vertexBufferDesc.CPUAccessFlags = 0;
	vertexBufferDesc.MiscFlags = 0;
	vertexBufferDesc.StructureByteStride = 0;

	// Give the subresource structure a pointer to the vertex data.
	vertexData.pSysMem = vertices;
	vertexData.SysMemPitch = 0;
	vertexData.SysMemSlicePitch = 0;

	// Now create the vertex buffer.
	result = device->CreateBuffer(&vertexBufferDesc, &vertexData, &m_vertexBuffer);
	if (FAILED(result))
	{
		return false;
	}

	// Set up the description of the static index buffer.
	indexBufferDesc.Usage = D3D11_USAGE_DEFAULT;
	indexBufferDesc.ByteWidth = sizeof(unsigned long) * m_indexCount;
	indexBufferDesc.BindFlags = D3D11_BIND_INDEX_BUFFER;
	indexBufferDesc.CPUAccessFlags = 0;
	indexBufferDesc.MiscFlags = 0;
	indexBufferDesc.StructureByteStride = 0;

	// Give the subresource structure a pointer to the index data.
	indexData.pSysMem = indices;
	indexData.SysMemPitch = 0;
	indexData.SysMemSlicePitch = 0;

	// Create the index buffer.
	result = device->CreateBuffer(&indexBufferDesc, &indexData, &m_indexBuffer);
	if (FAILED(result))
	{
		return false;
	}

	// Release the arrays now that the vertex and index buffers have been created and loaded.
	delete[] vertices;
	vertices = 0;

	delete[] indices;
	indices = 0;

	return true;
}

void Terrain::RenderBuffers(ID3D11DeviceContext * deviceContext)
{
	unsigned int stride;
	unsigned int offset;

	// Set vertex buffer stride and offset.
	stride = sizeof(VertexType);
	offset = 0;

	// Set the vertex buffer to active in the input assembler so it can be rendered.
	deviceContext->IASetVertexBuffers(0, 1, &m_vertexBuffer, &stride, &offset);

	// Set the index buffer to active in the input assembler so it can be rendered.
	deviceContext->IASetIndexBuffer(m_indexBuffer, DXGI_FORMAT_R32_UINT, 0);

	// Set the type of primitive that should be rendered from this vertex buffer, in this case triangles.
	deviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	return;
}
//Midpoint Displacement
float Terrain::jitter(float value, float spread)
{
	float r = (((float)rand() / (RAND_MAX)) * 2) - 1;
	return value + spread + (r * m_amplitude);
}

/*
  The four parameters define the square we're working on in the heightmap:
  lx: the x coordinate for the left-hand corners
  rx: the x coordinate for the right-hand corners
  by: the y coordinate for the bottom corners
  ty: the y coordinate for the top corners
*/
void Terrain::MidPointDisplacement(int lx, int rx, int by, int ty, float spread)
{
	int cx = (lx + rx) / 2;
	int cy = (by + ty) / 2;

	float topHeight = Average(m_heightMap[m_terrainHeight * ty + lx].y, m_heightMap[m_terrainHeight * ty + rx].y);
	int topIndex = m_terrainHeight * ty + cx;
	m_heightMap[topIndex].y = jitter(topHeight, spread);

	float leftHeight = Average(m_heightMap[m_terrainHeight * by + lx].y, m_heightMap[m_terrainHeight * ty + lx].y);
	int leftIndex = m_terrainHeight * cy + lx;
	m_heightMap[leftIndex].y = jitter(leftHeight, spread);

	float bottomHeight = Average(m_heightMap[m_terrainHeight * by + lx].y, m_heightMap[m_terrainHeight * by + rx].y);
	int bottomIndex = m_terrainHeight * by + cx;
	m_heightMap[bottomIndex].y = jitter(bottomHeight, spread);

	float rightHeight = Average(m_heightMap[m_terrainHeight * by + rx].y, m_heightMap[m_terrainHeight * ty + rx].y);
	int rightIndex = m_terrainHeight * cy + rx;
	m_heightMap[rightIndex].y = jitter(rightHeight, spread);

	float centreHeight = Average(topHeight, leftHeight, bottomHeight, rightHeight);
	int centreIndex = m_terrainHeight * cy + cx;
	m_heightMap[centreIndex].y = jitter(centreHeight, spread);
}
bool Terrain::GenerateHeightMap(ID3D11Device* device)
{
	bool result;

	//Initialise corners
	m_heightMap[0].y = ((float)rand() / (10000)) * m_amplitude;
	m_heightMap[m_terrainHeight * (m_terrainHeight - 1)].y = ((float)rand() / (10000)) * m_amplitude;
	m_heightMap[m_terrainWidth - 1].y = ((float)rand() / (10000)) * m_amplitude;
	m_heightMap[m_terrainHeight * (m_terrainHeight - 1) + m_terrainWidth - 1].y = ((float)rand() / (10000)) * m_amplitude;

	float spread = 0.3f;
	int exponent = log2(m_terrainWidth - 1); //Exponent will be 7 since heightmap's size is 129 (2^n + 1)
	//We need therefore 7 passes in order to complete our heightmap
	for (int i = 0; i < exponent; i++)
	{
		//At iteration i, we need to displace a 2^ix2^i collection of squares
		int chunks = pow(2, i);
		int chunkWidth = (m_terrainWidth - 1) / chunks;

		for (int x = 0; x < chunks; x++)
		{
			for (int y = 0; y < chunks; y++)
			{
				int leftX = chunkWidth * x;
				int rightX = leftX + chunkWidth;
				int bottomY = chunkWidth * y;
				int topY = bottomY + chunkWidth;
				MidPointDisplacement(leftX, rightX, bottomY, topY, spread);
			}
		}
		spread *= 0.5f;
		Smoothing(device);
		//Voronoi(device, 5);
	}

	result = CalculateNormals(); //Re-evaluate the normals after altering the vertices coordinates
	if (!result)
	{
		return false;
	}

	result = InitializeBuffers(device); //Re-initialize the buffers after altering the vertices coordinates
	if (!result)
	{
		return false;
	}
}

float Terrain::Average(float x, float y)
{
	return ((x + y) / 2.0f);
}

float Terrain::Average(float x, float y, float z, float w)
{
	return ((x + y + z + w) / 4.0f);
}


bool Terrain::Smoothing(ID3D11Device* device)
{
	int index;
	bool result;
	int neighbours[8];

	for (int j = 0; j < m_terrainHeight; j++)
	{
		for (int i = 0; i < m_terrainWidth; i++)
		{
			index = (m_terrainHeight * j) + i;

			//Top neighbours
			neighbours[0] = (m_terrainHeight * (j + 1)) + (i - 1);
			neighbours[1] = (m_terrainHeight * (j + 1)) + (i);
			neighbours[2] = (m_terrainHeight * (j + 1)) + (i + 1);

			//Side neighbours
			neighbours[3] = (m_terrainHeight * (j)) + (i - 1);
			neighbours[4] = (m_terrainHeight * (j)) + (i + 1);

			//Bottom neighbours
			neighbours[5] = (m_terrainHeight * (j - 1)) + (i - 1);
			neighbours[6] = (m_terrainHeight * (j - 1)) + (i);
			neighbours[7] = (m_terrainHeight * (j - 1)) + (i + 1);

			float sum = m_heightMap[index].y;
			for (int i = 0; i < 8; i++)
			{
				if (neighbours[i] >= 0 && neighbours[i] < m_terrainHeight * m_terrainWidth)
					sum += m_heightMap[neighbours[i]].y;
				else
					sum += m_heightMap[index].y; //if the index does not exist, replicate the current height value of the vertex
			}

			m_heightMap[index].y = sum / 9.0f;
		}
	}

	result = CalculateNormals(); //Re-evaluate the normals after altering the vertices coordinates
	if (!result)
	{
		return false;
	}

	result = InitializeBuffers(device); //Re-initialize the buffers after altering the vertices coordinates
	if (!result)
	{
		return false;
	}
}

bool Terrain::Voronoi(ID3D11Device* device, int nRegions)
{
	bool result;
	int* seedIndex = new int[nRegions];
	float height = 5.0f;

	for (int i = 0; i < nRegions; i++)
	{
		//Generate a random seed point
		int x = rand() % m_terrainWidth;
		int z = rand() % m_terrainHeight;
		seedIndex[i] = (z * m_terrainHeight) + x;
		//Check that the newly generated point does not exist already
		//...

		//Set the starting height values for the seed points
		m_heightMap[seedIndex[i]].y = height * i;
	}

	for (int j = 0; j < m_terrainHeight; j++)
	{
		for (int i = 0; i < m_terrainWidth; i++)
		{
			int index = (j * m_terrainHeight) + i;
			//Check that we're not on a seed point
			bool skipIteration = false;
			for (int k = 0; k < nRegions; k++)
			{
				if (index == seedIndex[k])
				{
					skipIteration = true;
					break;
				}
			}
			if (!skipIteration)
			{
				float min = m_terrainWidth * sqrt(2.0f); //Works only with square terrains
				for (int k = 0; k < nRegions; k++)
				{
					if (Distance(m_heightMap[index].x, m_heightMap[index].z, m_heightMap[seedIndex[k]].x, m_heightMap[seedIndex[k]].z) < min)
					{
						min = Distance(m_heightMap[index].x, m_heightMap[index].z, m_heightMap[seedIndex[k]].x, m_heightMap[seedIndex[k]].z);
						m_heightMap[index].y = m_heightMap[seedIndex[k]].y;
					}
				}
			}
		}
	}

	result = CalculateNormals(); //Re-evaluate the normals after altering the vertices coordinates
	if (!result)
	{
		return false;
	}

	result = InitializeBuffers(device); //Re-initialize the buffers after altering the vertices coordinates
	if (!result)
	{
		return false;
	}
}

bool Terrain::Update()
{
	return true; 
}

float Terrain::Distance(float lx, float lz, float rx, float rz)
{
	return sqrt(pow((rx - lx), 2) + pow((rz - lz), 2));
}
float Terrain::randomNumGen()
{
	float randomHeight = rand() % 10;
	return randomHeight;
}

float* Terrain::GetWavelength()
{
	return &m_wavelength;
}

float* Terrain::GetAmplitude()
{
	return &m_amplitude;
}
