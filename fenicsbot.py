import twitter
from time import sleep, time
from parser import parse, excise

class FEniCSbot(object):
    def __init__(self, twitter_api_object,
                 sleep_time=10, break_loop=lambda *_: False):
        self.api = twitter_api_object
        self.sleep_time = sleep_time
        self._last_check_id = 1
        self.break_loop = break_loop

    def print_status(self, message):
        """
        Helper function for printing a short status to standard out.

        :param message: Message to print.
        """

        print message

    def tweet_image(self, img_fn, tweet):
        """
        Tweets an image as a reply tweet.

        :param img_fn: Filename of image to tweet.
        :param tweet: Twitter status object of tweet to reply to.
        """

        print "-----------Tweeting a solution! ------------"
        solution_tweet_text = "@{}: I solved your problem ({})!".format(tweet.user.screen_name, excise(tweet.text).strip())[:100]
        self.api.PostMedia(solution_tweet_text, img_fn,
                           in_reply_to_status_id=str(tweet.id))

    def tweet_error(self, tweet, e):
        """
        Tweets an error message in reply to a tweet which did not parse.

        :param tweet: Tweet to reply to.
        :param e: The Exception object
        """
        self.print_status(".......FEniCSbot could not come to the rescue......")
        err_desc = "{}".format(e.message)

        doc_desc = " See docs: http://bit.ly/1T4KuNt"

        error_tweet = "@{}: I failed to solve your problem... ".format(tweet.user.screen_name) + err_desc + doc_desc

        print error_tweet
        self.api.PostUpdate(error_tweet[:140], in_reply_to_status_id=str(tweet.id))

    def get_mentions(self):
        """
        Gets tweets @FEniCSbot, and returns those which are recent, have not
        been answered yet, are not from @FEniCSbot, and include the word 'solve'.
        """

        new_mentions = self.api.GetSearch(term="@fenicsbot",
                                          since_id=self._last_check_id)
        if self._last_check_id == 1:
                # in first iteration, don't grab tweets from beginning of time
                is_recent = lambda t: (t.created_at_in_seconds >
                                       self.program_start_time)
                new_mentions = filter(lambda t: is_recent(t), new_mentions)

        if len(new_mentions) > 0:
            self._last_check_id = new_mentions[0].id


        # list of predicates which are False if a tweet should be ignored
        requirement_list = [
            # ignore tweets from self
            lambda t: t.user.screen_name.lower() != "fenicsbot",
            
            # FEniCSbot should not reply to tweets not including "solve"
            lambda t: "solve" in t.text.lower(),
            
            # FEniCSbot should not reply to retweets
            lambda t: hasattr(t, "retweeted_status")
        ]

        
        for check in check_list:
            new_mentions = filter(check, new_mentions)

        return new_mentions

    def start(self):
        """
        Starts FEniCSbot listen loop, causing it to
        reply to all tweets containing @FEniCSbot.
        """
        self.program_start_time = time()

        while not self.break_loop(self):
            self.print_status("\nScanning for new tweets...")

            for tweet in self.get_mentions():
                try:
                    solver = parse(tweet.text)
                    solver.solve()
                    img_fn = solver.plot()
                    self.tweet_image(img_fn, tweet)

                except Exception as e :
                    print "Got error ", e
                    print "Traceback: ",
                    import traceback
                    traceback.print_exc()
                    self.tweet_error(tweet, e)
                    self.print_status("Error replying to: {}".format(tweet.text))

            self.print_status("Scan for new tweets complete.\n")
            sleep(self.sleep_time)




if __name__ == "__main__":
    from secrets import secret_dict

    # twitter_api = twitter.Api(consumer_key=secret_dict["consumer_key"],
    #                           consumer_secret=secret_dict["consumer_secret"],
    #                           access_token_key=secret_dict["access_token_key"],
    #                           access_token_secret=secret_dict["access_token_secret"]
    #         )
    twitter_api = twitter.Api(**secret_dict) # should be equivalent to the above
    bot = FEniCSbot(twitter_api)

    bot.start()
