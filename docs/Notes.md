### **How to output point data**
#### **An example of 2D Point data format**
This is the document discussing how to output the simulation data in a general format so that the visualizer could easily visualize the results.

Please note that the raw data (every point) includes 4 digits:
```cpp
float data[4] = {position_x, position_y, reserved_1, reserved 2};
```
So the total memory would require:
```cpp
float *data = new float[number_of_points * 4];
```
**Note:** Of course you can design your own style to write out the points onto disk, following just an example of the idea to save current frame point locations. 
```cpp
void write_file(double *pos_x, double *pos_y, int num, int frame)
{
	char filename[256];
	sprintf(filename,"YourPath/Particle_data%04d.bin",frame);
	float *data;
	data = new float[num*4];
	for (int i=0;i<num;i++)
	{
		data[i*4+0] = pos_x[i];
		data[i*4+1] = pos_y[i];
		data[i*4+2] = 0;
		data[i*4+3] = 1;
	}
	FILE *f;
	f = fopen(filename,"wb");
	fwrite(&num, sizeof(int), 1, f);
	fwrite(data,sizeof(float),4*num,f);

	delete[]data;
	fclose(f);
}
```
