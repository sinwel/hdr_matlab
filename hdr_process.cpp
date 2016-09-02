

unsigned short	cacheline[12][4096];
unsigned short  dstline[8][4096];

void dma_cachetransf(unsigned short *frame, unsigned short *cache, int yoffset, int lines, int W, int stride_o, int stride_i)
{
	frame = frame + stride_o*yoffset;
	for (int y = 0; y < lines; y++)
	{
		for (int x = 0; x < W; x++)
		{
			cache[x] = frame[x];
		}
		cache += stride_i;
		frame += (stride_o);
	}
}

void dma_outtransf(unsigned short *frame, unsigned short *cache, int yoffset, int lines, int W, int stride_o, int stride_i)
{
	frame = frame + stride_o*yoffset;
	for (int y = 0; y < lines; y++)
	{
		for (int x = 0; x < W; x++)
		{
			frame[x] = cache[x];
		}
		cache += stride_i;
		frame += (stride_o);
	}
}

void hdr_block_process(unsigned short **pi, unsigned short **pd, int offset_pi, int offset_pd)
{
	for (int y = 0; y < 4; y++)
	{
		for (int x = 0; x < 32; x++)
		{
			pd[y][offset_pd + x] = pi[y + 2][offset_pi + x];
		}
	}
}

void hdrprocess_sony_raw(unsigned short *src, unsigned short *dst, unsigned short *simage, int W, int H, int times, int noise_thred)
{
	unsigned short	*psrc[12];
	unsigned short	*pdst[8];
	int				y_offset;

	for (int i = 0; i < 12; i++)
		psrc[i] = cacheline[i];

	for (int i = 0; i < 8; i++)
		pdst[i] = dstline[i];


	dma_cachetransf(src, psrc[2]+32, 0, 6, W, W, 4096);
	y_offset = 6;
	for (int y = 0; y < H; y += 4)
	{
		int	read_lines = 4;
		if (y_offset > H - 4)
		{
			read_lines = H - y_offset;
		}
		if (read_lines > 0)
		{
			dma_cachetransf(src, psrc[8]+32, y_offset, read_lines, W, W, 4096);
			y_offset += read_lines;
		}

		for (int x = 0; x < W; x += 32)
		{
			hdr_block_process(psrc, pdst, x+32, x);
		}

		dma_outtransf(dst, pdst[0], y, 4, W, W , 4096);


		for (int i = 0; i < 4; i++)
		{
			unsigned short *ptmp;

			ptmp = psrc[i];
			psrc[i] = psrc[i + 4];
			psrc[i + 4] = psrc[i + 8];
			psrc[i + 8] = ptmp;
		}

		for (int i = 0; i < 4; i++)
		{
			unsigned short *ptmp;

			ptmp = pdst[i];
			pdst[i] = pdst[i + 4];
			pdst[i + 4] = ptmp;
		}
	}




}