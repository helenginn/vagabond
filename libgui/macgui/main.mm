//
//  main.m
//  VagabondViewer
//
//  Created by Helen Ginn on 02/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "OptionSpawner.h"

int main(int argc, const char * argv[])
{
	OptionSpawner *spawner = [[OptionSpawner alloc] init];
	NSMutableArray *argvs = [NSMutableArray array];

	for (int i = 0; i < argc; i++)
	{
		[argvs addObject:[NSString stringWithCString:argv[i] encoding:NSUTF8StringEncoding]];
	}

	NSThread *thread = [[NSThread alloc] initWithTarget:spawner
											  selector:@selector(spawnVagabond:)
												object:argvs];
	[thread start];

	return NSApplicationMain(argc, argv);
}
