#pragma once

using namespace DirectX;

class Terrain
{
private:
	struct VertexType
	{
		DirectX::SimpleMath::Vector3 position;
		DirectX::SimpleMath::Vector2 texture;
		DirectX::SimpleMath::Vector3 normal;
	};
	struct HeightMapType
	{
		// hold the height map data. Each point on the height map will have an x, y, and z coordinate.
		float x, y, z;
		float nx, ny, nz;
		float u, v;
	};
public:
	Terrain();
	~Terrain();

	// transfers the array information into vertex buffers so that DirectX can render it.
	bool Initialize(ID3D11Device*, int terrainWidth, int terrainHeight);
	void Render(ID3D11DeviceContext*);
	bool GenerateHeightMap(ID3D11Device*);
	bool Update();
	float randomNumGen();
	float* GetWavelength();
	bool Smoothing(ID3D11Device*);
	float* GetAmplitude();
	bool Voronoi(ID3D11Device*, int nRegions);

private:
	bool CalculateNormals();
	void Shutdown();
	void ShutdownBuffers();
	float jitter(float, float);
	float Distance(float, float, float, float);

	void MidPointDisplacement(int, int, int, int, float);
	float Average(float, float);
	float Average(float, float, float, float);
	bool InitializeBuffers(ID3D11Device*);
	void RenderBuffers(ID3D11DeviceContext*);
	

private:
	bool m_terrainGeneratedToggle;
	int m_terrainWidth, m_terrainHeight;
	ID3D11Buffer * m_vertexBuffer, *m_indexBuffer;
	int m_vertexCount, m_indexCount;
	float m_frequency, m_amplitude, m_wavelength;
	HeightMapType* m_heightMap;

	//arrays for our generated objects Made by directX
	std::vector<VertexPositionNormalTexture> preFabVertices;
	std::vector<uint16_t> preFabIndices;
};

