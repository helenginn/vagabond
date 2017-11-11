//
//  OptionSpawner.m
//  VagabondViewer
//
//  Created by Helen Ginn on 03/11/2017.
//  Copyright © 2017 Strubi. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "OptionSpawner.h"
#include "Options.h"

@implementation OptionSpawner

- (void)spawnVagabond:(NSArray *)argvs
{
	char **argv = (char **)malloc(sizeof(char *) * [argvs count]);

	for (int i = 0; i < [argvs count]; i++)
	{
		long stringSize = [[argvs objectAtIndex:i] length] + 1;
		const char *oneArgv = [[argvs objectAtIndex:i] cStringUsingEncoding:NSUTF8StringEncoding];
		argv[i] = (char *)malloc(sizeof(char) * stringSize);
		memcpy(argv[i], oneArgv, stringSize * sizeof(char));
	}

	int count = (int)[argvs count];
	OptionsPtr myOptions = OptionsPtr(new Options(count, (const char **)argv));
	Options::setRuntimeOptions(myOptions);
	myOptions->run();
}

@end
