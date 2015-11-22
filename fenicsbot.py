import twitterp
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
        solution_tweet_text = "@{}: I solved your problem ({})!".format(tweet.user.screen_name, 
                                                                        excise(tweet.text).strip())[:100]
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

        error_tweet = "@{}: I failed to solve your problem... {} {}"
        error_tweet = error_tweet.format(tweet.user.screen_name, 
                                         err_desc, doc_desc)

        print error_tweet
        self.api.PostUpdate(error_tweet[:140], 
                            in_reply_to_status_id=str(tweet.id))
    def tweet_welcome(self, tweet):
        """
        Replies to a thank you tweet with something nice.
        
        :param tweet: Tweet to reply to.
        """
        welcome_message = "@{}: You're welcome!".format(tweet.user.screen_name)
        self.api.PostUpdate(welcome_messsage[:140], 
                            in_reply_to_status_id=str(tweet.id))
        
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
        check_list = [
            # ignore tweets from self
            lambda t: t.user.screen_name.lower() != "fenicsbot",
            
            # FEniCSbot should not reply to tweets not including "solve"
            lambda t: "solve" in t.text.lower(),
            
            # # FEniCSbot should not reply to retweets
            # lambda t: not hasattr(t, "retweeted_status")
            lambda t: t.text[:3] != "RT "
        ]
        
        solve_requests = new_mentions
        for check in check_list:
            solve_requests = filter(check, solve_requests)

        non_requests = filter(lambda t: t not in solve_requests, new_mentions)
        thank_yous = filter(lambda t: "thank" in t.text.lower(), non_requests)
        
        return solve_requests, thank_yous

    def start(self):
        """
        Starts FEniCSbot listen loop, causing it to
        reply to all tweets containing @FEniCSbot.
        """
        self.program_start_time = time()

        while not self.break_loop(self):
            self.print_status("\nScanning for new tweets...")
            
            solve_requests, thank_yous = self.get_mentions()
            for tweet in solve_requests:
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
            
            for tweet in thank_yous:
                self.tweet_welcome(tweet)

            self.print_status("Scan for new tweets complete.\n")
            sleep(self.sleep_time)




if __name__ == "__main__":
    from secrets import secret_dict

    twitter_api = twitter.Api(**secret_dict)
    bot = FEniCSbot(twitter_api)

    bot.start()
