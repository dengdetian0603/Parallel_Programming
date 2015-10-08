class shiftfill
{
	public void shiftfill()
	{
	}

	public static void main(String[] args)
	{
		long maxkey;

		for (int i=20; i<57; i++)
		{
			maxkey = ~(0L);
			maxkey = maxkey >>> (64-i);
			System.out.println("keysize bit: " + i + " keylength: " + maxkey);
		}
	}
}